from io import StringIO
from time import time
from sh import mafft
import logging

logger = logging.getLogger()


class MAFFTWrapper(object):  # global aligners, all vs all for MSA
    def __init__(self):
        pass

    def auto(
        self,
        fasta_path: str,
        output_file: str,
        gap_penalty: int = 4,
        extension_penalty: int = 1,
        threads: int = 30,
        buffer: bool = False,
        keep_length: bool = False
    ) -> StringIO:
        """Runs an alignment of all records in a fasta file.

        :param fasta_path: path to fasta
        :type fasta_path: str
        :param output_file: path to output
        :type output_file: str
        :param gap_penalty: gap penalty, defaults to 4
        :type gap_penalty: int, optional
        :param extension_penalty: extension penalty, defaults to 1
        :type extension_penalty: int, optional
        :param threads: threads, defaults to 30
        :type threads: int, optional
        :param buffer: skip file writting and return buffer, defaults to False
        :type buffer: bool, optional
        :return: buffered alignment or None
        :rtype: StringIO
        """
        if keep_length is True:
            keep_length = "--keeplength"
        else:
            keep_length = None

        tic = time()
        if buffer is False:
            with open(output_file, "w") as ofile:
                mafft(
                    "--auto",
                    "--op",
                    str(gap_penalty),
                    "--ep",
                    str(extension_penalty),
                    "--adjustdirection",
                    keep_length,
                    "--preservecase",
                    "--treeout",
                    "--thread",
                    str(threads),
                    fasta_path,
                    _out=ofile,
                )
        else:
            buf = StringIO()
            mafft(
                "--auto",
                "--op",
                str(gap_penalty),
                "--ep",
                str(extension_penalty),
                " --adjustdirection",
                keep_length,
                "--preservecase",
                "--treeout",
                "--thread",
                str(threads),
                fasta_path,
                _out=buf,
            )
            buf.seek(0)
        elapsed_time = time() - tic
        message = f"Time elapsed in mafft alignment: {elapsed_time}"

        logger.info(message)
        if buffer is True:
            return buf

    def auto_add_fragments(
        self,
        fasta_path: str,
        base_fasta: str,
        output_file: str,
        gap_penalty: int = 4,
        extension_penalty: int = 1,
        threads: int = 30,
        buffer: bool = False,
        keep_length: bool = True

    ) -> StringIO:
        """Automatic alignment against a reference

        :param fasta_path: path to fasta
        :type fasta_path: str
        :param ref_file: path to reference fasta
        :type ref_file: str
        :param output_file: output file
        :type output_file: str
        :param gap_penalty: gap penalty, defaults to 4
        :type gap_penalty: int, optional
        :param extension_penalty: extension penalty, defaults to 1
        :type extension_penalty: int, optional
        :param threads: threads for analysis, defaults to 30
        :type threads: int, optional
        :param buffer: skip file writting and return buffer, defaults to False
        :type buffer: bool, optional
        :return: aligned buffer or None
        :rtype: StringIO
        """
        if keep_length is True:
            keep_length = "--keeplength"
        else:
            keep_length = None
        tic = time()
        if buffer is False:
            with open(output_file, "w") as ofile:
                mafft(
                    "--auto",
                    "--op",
                    str(gap_penalty),
                    "--ep",
                    str(extension_penalty),
                    "--adjustdirection",
                    "--preservecase",
                    keep_length,
                    "--thread",
                    str(threads),
                    "--addfragments",
                    fasta_path,
                    base_fasta,
                    _out=ofile,
                )
        else:
            buf = StringIO()
            mafft(
                "--auto",
                "--op",
                str(gap_penalty),
                "--ep",
                str(extension_penalty),
                "--adjustdirection",
                "--preservecase",
                keep_length,
                "--thread",
                str(threads),
                "--addfragments",
                fasta_path,
                base_fasta,
                _out=buf,
            )
        elapsed_time = time() - tic
        message = f"Time elapsed in mafft alignment: {elapsed_time}"

        logger.info(message)
        if buffer is True:
            return buf