import logging
import os
import re
import warnings
from pathlib import Path

import CIAlign.utilityFunctions as cuf
from black import out
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord

from .plot import draw_mini_alignment, plot_pssm
from .wrappers import MAFFTWrapper
#from .parser import FastaParse

warnings.filterwarnings("ignore")
mafftwrapper = MAFFTWrapper()

logger = logging.getLogger(__name__)

# This entire file is a complete mess and requires a hard refactor

class MSABase(object):
    """Multiple sequence alignment shared functions. Able to align seqs, generate pssm and guides.

    Args:
        object object): Base form.
    """

    def __init__(self, parser: FastaParse):
        self.parser = parser
        self.msa_fasta_path = None
        self.expanded_msa_fasta_path = None
        self.segement_length = None
        self.pssm_df = None
        self.negative_taxid_path = None
        self.local_fda_path = None
        self.input_sequences_path = None

    def base_file_name(self):
        return self.parser.base_file_name

    def get_preset(self):
        return self.parser.preset


    def get_pssm_perc(self):
        preset = self.get_preset()
        return float(preset.pssm_perc)

    def get_curated_fasta(self):
        return self.parser.curated_fasta_path

    def get_taxids(self):
        return self.parser.preset.negative_taxids

    def get_uncurated_fasta(self):
        return self.parser.fasta_path

    def get_subject_fasta(self):
        outlier_fasta = self.parser.outlier_curated_fasta_path
        if outlier_fasta is not None:
            return outlier_fasta
        else:
            return self.parser.fasta_path

    def get_representative_path(self):
        return self.parser.representative_fasta_path

    def set_representative_path(self, new_ref):
        """Set project representative path

        Args:
            new_ref (str]): path to representative
        """
        if os.path.exists(new_ref):
            self.parser.representative_fasta_path = new_ref
        else:
            message = f"File: {new_ref} does not exist."

            logger.warning(message)

    def get_filtered_refs_df(self):
        return self.parser.filtered_refs

    def get_preset_adjusted(self):
        preset = self.get_preset()
        return float(preset.blast_adjusted_perc)

    def get_params(self):
        """Function to quickly retrieve paramaters df

        Returns:
            pandas.DataFrame: dataframe containing paramaters
        """
        dic = self.parser.preset.__dict__
        dic.update(self.parser.__dict__)
        dic.update(self.__dict__)
        bad_terms = [
            "input_refs",
            "filtered_refs",
            "branch_df",
            "negative_taxid_path",
            "alignment_array",
            "alignment_ids",
            "pssm_df",
            "parser",
            "preset",
            "consensus_sequence",
            "dumb_consensus",
            "gapless_df",
            "gap_consensus",
            "features_df",
        ]
        dic["preset_type"] = str(self.get_preset())
        for x in bad_terms:
            if x in dic:
                dic.pop(x)
        dic["negative_taxids"] = ",".join([str(tax) for tax in dic["negative_taxids"]])
        for cat in dic:
            val = dic[cat]
            if not isinstance(val, str):
                dic[cat] = str(dic[cat])
        df = pd.DataFrame(dic, index=["input"]).T.dropna()
        return df

    def save_params(self):
        """Save parameters to local file."""
        df = self.get_params()
        df.to_csv(f"{self.base_file_name()}_parameters.csv")

    def get_accession_id(self):
        """
        Pull accession id from fasta header of representative file

        Returns:
            str: accession id
        """
        rep_path = self.get_representative_path()
        acc_id = next(SeqIO.parse(rep_path, "fasta")).id
        return acc_id

    def extract_features_from_alignment(
        self,
        feature_name: str,
        feature_seq: str = None,
        feature_size_cutoff: float = 0.75,
    ) -> tuple:
        """extract a feature from an MSA

        :param feature_name: name of feature, should match gbff feature name
        :type feature_name: str
        :param feature_seq: custom feature seq to look for, defaults to None
        :type feature_seq: str, optional
        :param feature_size_cutoff: ungapped feature seq must be within certain size, defaults to 0.75
        :type feature_size_cutoff: float, optional
        :return: filepaths to (ungapped fasta, gapped/aligned fasta)
        :rtype: tuple
        """

        feature_size_cutoff_upper = 1 + (1 - feature_size_cutoff)
        aln_arr, aln_ids = self.alignment_array, self.alignment_ids
        rep_seq = "".join(aln_arr[0])

        feat_seq_pat = "-*".join(feature_seq)
        r = re.search(feat_seq_pat, rep_seq)
        if not r:
            feat_seq_pat = "-*".join(Seq(feature_seq).reverse_complement().__str__())

            r = re.search(feat_seq_pat, rep_seq)
        if r:
            start, end = r.span()

        sub_arr = aln_arr[:, start:end]

        recs = []
        aln_recs = []
        feat_len = len(feature_seq)
        passing_records = 0
        for i, row in zip(aln_ids, sub_arr):
            row_seq = "".join(row)
            row_seq_ungapped = re.sub("-*", "", row_seq)
            row_seq_size = len(row_seq_ungapped)
            row_seq_perc = row_seq_size / feat_len
            if (
                row_seq_perc >= feature_size_cutoff
                and row_seq_perc <= feature_size_cutoff_upper
            ):
                srec_seq = Seq(row_seq_ungapped)
                srec = SeqRecord(
                    id=i.split(" ")[0], description=i + " " + feature_name, seq=srec_seq
                )
                recs.append(srec)

                aln_srec_seq = Seq(row_seq)
                aln_srec = SeqRecord(
                    id=i.split(" ")[0],
                    description=i + " " + feature_name,
                    seq=aln_srec_seq,
                )
                aln_recs.append(aln_srec)
                passing_records += 1

        passed_perc = round(passing_records / len(sub_arr) * 100, 3)
        logger.info(
            f"Records which have feature {feature_name} within feature_size_cutoff of {feature_size_cutoff}: {passed_perc}%"
        )

        feat_fasta = f"{self.parser.base_file_name}_{feature_name}_features.fasta"
        feat_fasta_aln = (
            f"{self.parser.base_file_name}_{feature_name}_features_MAFFT_self.fasta"
        )

        SeqIO.write(recs, feat_fasta, "fasta")
        SeqIO.write(aln_recs, feat_fasta_aln, "fasta")

        return feat_fasta, feat_fasta_aln

    @staticmethod
    def get_complement_for_row(row: list) -> list:
        """Function to assist in calculating the RC of an alignment.
        optimized for speed.

        :param row: Aligned sequence
        :type row: list
        :return: RC aligned sequennce
        :rtype: list
        """
        mapt = {
            "A": "T",
            "G": "C",
            "T": "A",
            "C": "G",
            ".": ".",
            "-": "-",
            "N": "N",
            "R": "Y",
            "Y": "R",
            "K": "M",
            "M": "K",
            "S": "S",
            "W": "W",
            "B": "V",
            "D": "H",
            "H": "D",
            "V": "B",
            "X": "X",
        }
        nrow = []
        for i, x in enumerate(row):
            x = x.upper()
            y = mapt[x]
            nrow.append(y)

        return nrow



class MafftMSA(MSABase):
    """
    Class for the handling of mafft MSA features including: running mafft for a MSA output,
    loading an existing MSA, scoring the consensus threshold
    """

    def __init__(self, parser: FastaParse = FastaParse()):
        super().__init__(parser)
        self.alignment_array = None
        self.reverse_alignment_array = None
        self.alignment_ids = None
        self.consensus_sequence = None
        self.gapless_df = None
        self.gap_consensus = None
        self.dumb_consensus = None
        self.pssm_df = None
        self.all_tile_scores = None
        self.degenerate_consensus = None

        logger.info("Mafft MSA loaded parser")

    def reverse_complement_msa_array(self, aln_arr: np.array) -> np.array:
        """Calculate the reverse complement from a numpy array.
        Optimized for speed.

        :param msa_df: aligned sequences
        :type aln_arr: np.array
        :return: reverse complement of aligned sequences
        :rtype: np.array
        """

        aln_arr = self.alignment_array
        rmsa = np.apply_along_axis(self.get_complement_for_row, 1, aln_arr)
        rmsa = rmsa[:, ::-1]

        return rmsa

    def _fragments_has_outliers(self, msa_fasta_path: str):
        """
        Checks the newick tree of a mafft alignment that was run WITH addfragments.
        The phylogeny must be calculated and may take time as long as it took to originally
        run the addfragments method.
        """
        phy = Phylogeny()
        tree_path, tree, branch_df = phy.build_tree_from_alignment(msa_fasta_path)
        self.branch_df = branch_df
        outliers = phy.get_outliers(branch_df)
        if outliers:
            message = f"Outliers found in {msa_fasta_path}: {outliers}. Recalculating MSA unless remove_outliers=False."

            logger.info(message)
            parser = self.parser
            parser.filter_header_term(outliers, plot=False)
            outlier_path = f"outlier_{self.get_curated_fasta()}"
            parser.outlier_curated_fasta_path = outlier_path
            os.rename(self.get_curated_fasta(), outlier_path)
            parser.write_curated_fasta()
            return True
        else:
            message = f"No outliers found in {msa_fasta_path}."

            logger.info(message)
            return False

    def add_fragments_to_existing(
        self,
        additional_fasta: str,
        existing_fasta: str = None,
        output_msa_path: str = None,
        load_alignment: bool = True,
        keep_length: bool = True,
        plot: bool = True,
    ):

        op = self.get_preset().mafft_op
        ep = self.get_preset().mafft_ep
        threads = self.parser.preset.threads

        if output_msa_path is None:
            output_msa_path = self.base_file_name() + "_MAFFT_MSA_with_fragments.fasta"
        if existing_fasta is None and self.msa_fasta_path is not None:
            existing_fasta = self.msa_fasta_path
        elif existing_fasta is None and self.msa_fasta_path is None:
            raise Exception(
                "Cannot align additonal_fasta to existing_fasta. File not found"
            )

        mafftwrapper.auto_add_fragments(
            additional_fasta,
            existing_fasta,
            output_msa_path,
            gap_penalty=op,
            extension_penalty=ep,
            threads=threads,
            keep_length=keep_length,
        )

        if load_alignment is True:
            self.set_alignment(
                msa_fasta=output_msa_path,
                plot=plot,
                title=self.base_file_name() + "_partial",
            )
        self.expanded_msa_fasta_path = output_msa_path

    def mafft_alignment(
        self,
        fasta_path: str = None,
        rep: str = True,
        remove_outliers: bool = None,
        threads: int = None,
        keep_length: bool = False,
    ):
        """TODO This needs to be cleaned up."""
        if fasta_path is None:
            fasta_path = self.get_curated_fasta()
        stem = Path(fasta_path).stem
        self.msa_fasta_path = f"{stem}_MAFFT_MSA.fasta"
        if threads is None:
            threads = self.parser.preset.threads
        if remove_outliers is None:
            remove_outliers = self.get_preset().remove_outliers
        # If no representative file is set, then the standa rd mafft alignment will be computed.
        if rep is None:
            op = self.get_preset().mafft_op
            ep = self.get_preset().mafft_ep
            mafftwrapper.auto(
                fasta_path,
                self.msa_fasta_path,
                threads=threads,
                gap_penalty=op,
                extension_penalty=ep,
                keep_length=keep_length,
            )

            # if the guide tree exists, loading it is computationally cheap, but is not always
            # as accurate as the recalculated identity matrix
            # tree_path = f'{self.get_curated_fasta()}.tree'
            # self.tree_path = tree_path
            if remove_outliers is True:
                status = self._fragments_has_outliers(self.msa_fasta_path)
                if status is True:
                    mafftwrapper.auto(
                        fasta_path,
                        self.msa_fasta_path,
                        gap_penalty=op,
                        extension_penalty=ep,
                        threads=threads,
                        keep_length=keep_length,
                    )

        # If a representative is set, a decision tree will either run fragments with the representative
        # or choose the longest sequence as the representative to perform the alignment with.
        else:
            op = self.get_preset().mafft_op
            ep = self.get_preset().mafft_ep
            message = "MAFFT alignment using fragments setting and representative"
            logger.info(message)

            rep_path = self.get_representative_path()

            # Assumes a path was given
            if isinstance(rep, str):
                if os.path.exists(rep):
                    self.set_representative_path(rep)
                    rep_path = rep
                    mafftwrapper.auto_add_fragments(
                        fasta_path,
                        rep,
                        self.msa_fasta_path,
                        gap_penalty=op,
                        extension_penalty=ep,
                        threads=threads,
                        keep_length=keep_length,
                    )
                else:
                    message = f"Representative path {rep} not found, defaulting to longest seq."

                    logger.warning(message)
                    self.mafft_alignment(rep=True, threads=threads)

            # Assumes rep set to true, and a representitave was already set at some point in the pipeline
            elif rep_path is not None:
                # sets the representative path
                mafftwrapper.auto_add_fragments(
                    fasta_path,
                    rep_path,
                    self.msa_fasta_path,
                    gap_penalty=op,
                    extension_penalty=ep,
                    threads=threads,
                    keep_length=keep_length,
                )

            # Rewrite the longest sequence and write it to a file to set as the representative
            elif rep_path is None:
                parser = self.parser
                parser.write_curated_fasta(sep_ref=True)
                mafftwrapper.auto_add_fragments(
                    fasta_path,
                    rep_path,
                    self.msa_fasta_path,
                    gap_penalty=op,
                    extension_penalty=ep,
                    threads=threads,
                    keep_length=keep_length,
                )

            # The msa must be have a full phylogenetic distance identity matrix calculation done to determine
            # outliers, but only once. If there are outliers, the msa must be recomputed.
            if remove_outliers is True:
                status = self._fragments_has_outliers(self.msa_fasta_path)
                if status is True:
                    mafftwrapper.auto_add_fragments(
                        fasta_path,
                        rep_path,
                        self.msa_fasta_path,
                        gap_penalty=op,
                        extension_penalty=ep,
                        threads=threads,
                    )
        self.set_alignment()

    def auto(self):
        """
        Function to automatically run mafft alignment to guide sensitivity generation
        """
        self.mafft_alignment(rep=True)
        self.write_universal_regions()
        preset = self.get_preset()
        if preset.platform in ["ISR", "inspectr"] and "DNA" in preset.seq_type:
            self.guide_features_annotate()
        fasta_path = self.get_subject_fasta()
        self.blastutils.blast_sequences_sensitivity(
            self.input_sequences_path, fasta_path=fasta_path
        )

    def set_alignment(self, msa_fasta=None, plot=None, title=None):
        if msa_fasta is None:
            msa_fasta = self.msa_fasta_path
        else:
            self.msa_fasta_path = msa_fasta
        if title is not None:
            pass
        elif self.parser is not None:
            title = self.parser.base_file_name
        else:
            title = Path(msa_fasta).stem
        self.alignment_array, self.alignment_ids = cuf.FastaToArray(msa_fasta)
        if plot is None:
            plot = self.get_preset().plot_msa
        if plot is True:
            draw_mini_alignment(self.alignment_array, self.alignment_ids, title)

    def calculate_pssm(
        self, msa_fasta_path: str = None, pssm_perc: float = None, plot: bool = True
    ) -> tuple:
        if pssm_perc is None:
            pssm_perc = self.get_pssm_perc()
        else:
            self.get_preset().pssm_perc = pssm_perc
        if pssm_perc > 1.0:
            pssm_perc /= 100
        if msa_fasta_path is None:
            msa_fasta_path = self.msa_fasta_path
        c_align = next(AlignIO.parse(msa_fasta_path, "fasta"))
        summary_align = AlignInfo.SummaryInfo(c_align)
        rois, pssm_df, template, all_tile_scores = score_alignment_pssm(
            summary_align,
            self.alignment_array[0],
            pssm_perc=pssm_perc,
            tile_size=self.get_tile_window(),
        )
        self.pssm_df = pssm_df
        self.all_tile_scores = all_tile_scores
        tdf = pd.DataFrame(all_tile_scores)
        if plot is True:
            plot_pssm(tdf, self.parser.base_file_name)

        return rois, pssm_df, all_tile_scores

    def calculate_degenerate_sequence(
        self,
        pssm: pd.DataFrame = None,
        consensus_thresh: float = 95.0,
        require_gaps_below_consensus: bool = True,
    ) -> str:
        """automatically calculates the degenerate sequence from a
        MSA

        :param pssm: pssm dataframe, defaults to None
        :type pssm: pd.DataFrame, optional
        :param consensus_thresh: threshold needed to create degenerates for, defaults to 95.0
        :type consensus_thresh: float, optional
        :param require_gaps_below_consensus: , defaults to True
        :type require_gaps_below_consensus: bool, optional
        :return: degenerate consensus string
        :rtype: str
        """

        if pssm is not None:
            pass
        elif self.pssm_df is not None:
            pssm = self.pssm_df
        else:
            logger.warning("Missing PSSM calculation. Calculating PSSM.")
            self.calculate_pssm()
            pssm = self.pssm_df

        if consensus_thresh > 1:
            consensus_thresh /= 100

        non_gap_limit = 1 - consensus_thresh
        inds = pssm.columns.to_numpy()

        nuc = ["A", "T", "G", "C", "U"]

        degen_dict = {
            "AG": "R",
            "T": "T",
            "CT": "Y",
            "CG": "S",
            "AT": "W",
            "GT": "K",
            "AC": "M",
            "CGT": "B",
            "AGT": "D",
            "ACT": "H",
            "ACG": "V",
            "X": "X",
            "A": "A",
            "G": "G",
            "C": "C",
            "ACGT": "N",
        }

        msa_degenerate = ""

        for row in pssm.to_numpy():
            row = np.array([inds, row])
            row = row.T[np.argsort(row[1, :])]
            row = np.flip(row, axis=0)
            cum_sum = row[:, 1].cumsum()
            row = np.vstack((row.T, cum_sum))
            row_exceeds_gap_tolerance = False
            for i, (cu, letter, weight) in enumerate(zip(row[2], row[0], row[1])):
                if cu >= consensus_thresh:
                    break
                if require_gaps_below_consensus is True and letter == "-":
                    if weight > non_gap_limit:
                        msa_degenerate += "-"
                        row_exceeds_gap_tolerance = True
            if row_exceeds_gap_tolerance is True:
                continue
            row = row[:, : i + 1]
            if row[0][0] == "-":
                msa_degenerate += "-"
                continue
            nuc_key = np.array([nu for nu in row[0] if nu in nuc])
            nuc_key.sort()
            if any(nuc_key):
                nuc_key = (
                    nuc_key.astype(object)
                    .sum()
                    .upper()
                    .replace("U", "T")
                    .replace("-", "")
                )
                degen = degen_dict[nuc_key]
            else:
                degen = "-"

            msa_degenerate += degen

        self.degenerate_consensus = msa_degenerate

        return msa_degenerate

  