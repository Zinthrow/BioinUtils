class EntrezWrapper(object):
    """
    Downloads files from the NCBI nucleotide database.

    Input:
    email: str, entrez requires an email for every query
    ncbi_api_key: str, strongly reccomended to include this for large queries
    """

    def __init__(self, email, ncbi_api_key=None):
        self.entrez = Entrez
        self.entrez.email = email
        if ncbi_api_key is not None:
            self.entrez.api_key = ncbi_api_key
        self.entrez.max_tries = 10
        self.entrez.sleep_between_tries = 10

    def get_search_dict(
        self, database: str, mquery: str, sort: str = "significance", retmax: int = None
    ) -> tuple:
        """Retrieve the initial search id hit result.

        :param database: ncbi database path local (path/to/database/nt)
        :type database: str
        :param mquery: query term
        :type mquery: str
        :param sort: sortby, defaults to 'significance'
        :type sort: str, optional
        :param retmax: top max retry, defaults to None
        :type retmax: int, optional
        :return: search results, search count
        :rtype: tuple
        """
        handle = self.entrez.esearch(
            database, term=mquery, usehistory="y", sort=sort, retmax=retmax
        )
        search_results = Entrez.read(handle)
        s_count = int(search_results["Count"])
        message = (
            f"Found {s_count} entries associated with {mquery} in {database} database"
        )
        logger.info(message)

        return search_results, s_count

    def get_elink_dict(
        self, id_term: str, db: str = "nuccore", target: str = "assembly"
    ) -> tuple:
        """
        Retrieve the elink associated with a accession number to retrieve the corresponding
        assembly entry. Allows the program to use contig and scaffolding accession ids if
        no master accession id can be found.

        :param id_term: search term
        :type id_term: str
        :param db: database, defaults to 'nuccore'
        :type db: str, optional
        :param target: database target, defaults to 'assembly'
        :type target: str, optional
        :return: search results, search count
        :rtype: tuple
        """

        buffer = StringIO()
        elink(id=id_term, db=db, target=target, _out=buffer)
        buffer.seek(0)
        xml_string = buffer.read()
        search_results = xmltodict.parse(xml_string)["ENTREZ_DIRECT"]
        s_count = int(search_results["Count"])

        return search_results, s_count

    def parse_taxa_result(self, tdict: dict, query: str) -> dict:
        """retrieve and create taxonomic lineage

        :param tdict: search results output
        :type tdict: dict
        :param query: original query
        :type query: str
        :return: lineage
        :rtype: dict
        """

        if "eFetchResult" in tdict:
            message = f"No taxonomy information found concerning: {query}"

            logger.warning(message)
            return {"No Results": {"ScientificName": None, "TaxId": None}}

        tdict = tdict["TaxaSet"]["Taxon"]
        tname = tdict["ScientificName"]
        trank = tdict["Rank"]
        ttax = tdict["TaxId"]

        lineage = dict()
        for subdic in tdict["LineageEx"]["Taxon"]:
            name = subdic["ScientificName"]
            rank = subdic["Rank"]
            taxid = subdic["TaxId"]
            lineage[rank] = {"ScientificName": name, "TaxId": taxid}

        lineage[trank] = {"ScientificName": tname, "TaxId": ttax}

        return lineage

    def get_taxa_lineage(self, query: str) -> dict:
        """
        Retrieve the full taxonomic lineage

        :param query: query of organism
        :type query: str
        :return: organism taxonomic lineage
        :rtype: dict
        """
        database = "taxonomy"
        search_results, s_count = self.get_search_dict(database, query, retmax=1)
        if s_count == 0:
            search_results, s_count = self.get_elink_dict(query, target="taxonomy")
        fetch_handle = self.entrez.efetch(
            db="taxonomy",
            rettype="xml",
            retmode="xml",
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )
        fetch_xml = fetch_handle.read()
        tdict = xmltodict.parse(fetch_xml)
        lineage = self.parse_taxa_result(tdict, query)

        return lineage

    def get_nucleotide_fasta(
        self,
        query: str,
        outfile: str,
        q_type: str = "partial",
        db: str = "nucleotide",
        write: bool = True,
        retmax: int = 5000,
    ) -> str:
        """
        :param query: organism query
        :type query: str
        :param outfile: file to write output to
        :type outfile: str
        :param q_type: [genome, plasmid, custom], defaults to 'genome'
        :type q_type: str, optional
        :param db: database in ncbi, defaults to 'nucleotide'
        :type db: str, optional
        :param write: write to file, defaults to True
        :type write: bool, optional
        :param retmax: how many records to pull per query, defaults to 5000
        :type retmax: int, optional
        :raises ValueError: No results associate with query
        :return: results
        :rtype: str
        """
        query = query.replace("_", " ")
        if q_type == "genome":
            mquery = f"{query} complete genome AND {query}[Organism]"
        elif q_type == "partial":
            mquery = f"{query} AND {query}[Organism]"
        elif q_type == "plasmid":
            mquery = f"{query}[Title] AND plasmid[filter] AND {query}[Organism]"
        elif q_type == "custom":
            mquery = query
        else:
            mquery = f"{query}[Title]"

        search_results, s_count = self.get_search_dict(db, mquery)

        if s_count == 0:
            message = f"Query invalid: {query}"
            logger.error(message)
            raise ValueError(message)

        master_fetch = ""
        for start in range(0, s_count, retmax):
            fetch_handle = self.entrez.efetch(
                db=db,
                rettype="fasta",
                retmode="text",
                retstart=start,
                retmax=retmax,
                webenv=search_results["WebEnv"],
                query_key=search_results["QueryKey"],
            )
            master_fetch += fetch_handle.read()
            master_fetch += "\n"
        if write is True:
            with open(outfile, "w") as fasta_file:
                fasta_file.write(master_fetch)
        else:
            return master_fetch

    def format_assembly_params(
        self, query: str, assembly_db: str, out_dir: str, only_reps: bool, q_type: str
    ) -> tuple:
        """Add format for assembly database

        :param query: organism query
        :type query: str
        :param assembly_db: [genbank, refseq]
        :type assembly_db: str
        :param out_dir: output directory
        :type out_dir: str
        :param only_reps: representatives only
        :type only_reps: bool
        :param q_type: query type
        :type q_type: str
        :return: query, modified_query, output directory
        :rtype: tuple
        """
        query = query.replace(" ", "+").replace("_", "+")

        if q_type == "genome":
            query = query.replace("_", " ")
            mquery = (
                f"{query} AND latest+{assembly_db}[filter] AND complete+genome[filter]"
            )
        elif q_type == "partial":
            query = query.replace("_", " ")
            mquery = f'{query} AND latest+{assembly_db}[filter] AND ("complete genome"[filter] OR "chromosome level"[filter] OR "scaffold level"[filter])'
        else:
            mquery = query

        # First attempt to find a representative or reference genome otherwise pull
        if only_reps is True:
            qref = " AND (representative+genome[filter] OR reference+genome[filter])"
            mquery = mquery + qref
        else:
            mquery = mquery

        # Make out_dir if out_dir specified and does not exist and format
        if out_dir:
            if "/" not in out_dir:
                out_dir = out_dir + "/"
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

        return query, mquery, out_dir

    def write_merged_assembly(
        self, query: str, genome_list: list, bstr: bytes, only_reps: bool
    ) -> str:
        """
        Write assembly files merged into one fasta file.

        :param query: original query
        :type query: str
        :param genome_list: list of assembly genome files
        :type genome_list: list
        :param bstr: bytes fasta
        :type bstr: bytes
        :param only_reps: only use representatives
        :type only_reps: bool
        :return: fasta file name
        :rtype: str
        """

        if only_reps is True or len(genome_list) == 1:
            fname = genome_list[0]
        else:
            fname = query.replace("+", "_") + "_records.fasta"
        f = open(fname, "w")
        f.write(bstr)
        logger.info(f"Concatenated assembly fastas written at {fname}.")
        return fname

    def get_master_summary(
        self, s_count: int, retmax: int, search_results: dict
    ) -> list:
        """Get list of master summaries

        :param s_count: number of search results
        :type s_count: int
        :param retmax: records to pull at once
        :type retmax: int
        :param search_results: all search results
        :type search_results: dict
        :raises ValueError: invalid document summary.
        :return: master summary
        :rtype: list
        """
        master_summary = []
        for start in range(0, s_count, retmax):
            # Pull summary for all esearch results
            summary_handle = self.entrez.esummary(
                db="assembly",
                rettype="xml",
                retmode="xml",
                retstart=start,
                retmax=retmax,
                webenv=search_results["WebEnv"],
                query_key=search_results["QueryKey"],
            )

            summary = self.entrez.read(summary_handle, validate=False)
            # Load summary into dict
            if "DocumentSummarySet" in summary:
                summary = summary["DocumentSummarySet"]["DocumentSummary"]
                master_summary.extend(summary)
            else:
                logger.error(summary)
                raise ValueError(summary)

        return master_summary

    def get_assembly_fastas(
        self,
        query: str,
        out_dir: str = "",
        only_reps: bool = False,
        assembly_db: str = "genbank",
        q_type: str = "genome",
        read_only: bool = False,
        limit: int = None,
        individual: bool = True,
        merge_contigs: bool = False,
        retmax: int = 5000,
    ) -> list:
        """Function to retrieve all assembly fastas related to term.

        :param query: organism search query
        :type query: str
        :param out_dir: diretory to output fastas, defaults to ''
        :type out_dir: str, optional
        :param only_reps: representatives only, defaults to False
        :type only_reps: bool, optional
        :param assembly_db: database to pull assembly files from [genbank, refseq], defaults to 'genbank'
        :type assembly_db: str, optional
        :param q_type: query type [genome, partial], defaults to 'genome'
        :type q_type: str, optional
        :param read_only: , defaults to False
        :type read_only: bool, optional
        :param limit: only pull down to N records, defaults to None
        :type limit: int, optional
        :param individual: write individual assembly files, defaults to True
        :type individual: bool, optional
        :param merge_contigs: merge contigs into singular fasta header, defaults to False
        :type merge_contigs: bool, optional
        :param retmax: number of records to grab at once, defaults to 5000
        :type retmax: int, optional
        :raises ValueError: invalid query
        :return: list of genome records
        :rtype: list
        """

        # modify the query to suit given kwargs, critical step
        query, mquery, out_dir = self.format_assembly_params(
            query, assembly_db, out_dir, only_reps, q_type
        )

        database = "assembly"
        search_results, s_count = self.get_search_dict(database, mquery)
        if read_only is True:
            return search_results, s_count

        if s_count == 0:
            message = f"Query invalid: {mquery}"
            logger.error(message)

            raise ValueError(message)

        master_summary = self.get_master_summary(s_count, retmax, search_results)

        # Iterate through all results
        genome_list = []
        for s in master_summary[:limit]:
            path = s["FtpPath_GenBank"]
            stem = path.split("/")[-1]
            filepath = f"{stem}_genomic.fna.gz"
            fna_ftp = f"{path}/{filepath}"
            if out_dir:
                filepath = out_dir + filepath
            call(["wget", fna_ftp, "-O", filepath], stderr=subprocess.DEVNULL)
            call(["gunzip", filepath, "-f"])
            genome_list.append(filepath.rstrip(".gz"))

        if individual is False:
            bstr = ""
            for fil in genome_list:
                if not os.path.exists(fil):
                    continue
                file_string = open(fil, "r").read()
                if merge_contigs is True:
                    fil2 = re.sub(">.+", "\nXXXXXXXXXXXXXXXXXXXXXXX\n", file_string)
                    header = re.search(">.+", file_string).group()
                    file_string = header + fil2[3:]
                bstr += file_string
            fname = self.write_merged_assembly(query, genome_list, bstr, only_reps)
            return fname
        return genome_list

    def get_gbff_handle(self, acc_id: str, write: bool = False):
        """
        Get nucleotide gbff file

        Input:
        acc_id: str, ncbi accession id from fasta header.
        """
        database = "nucleotide"
        search_results, s_count = self.get_search_dict(database, acc_id)

        if s_count == 0:
            message = f"Query invalid: {acc_id}"
            logger.error(message)
            raise ValueError(message)

        fetch_handle = self.entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retstart=0,
            retmax=1,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )

        if write is False:
            return fetch_handle
        else:
            with open(f"{acc_id}.gbff", "w") as fasta_file:
                fasta_file.write(fetch_handle.read())

    def get_assembly_gbff(self, acc_id: str, limit: int = 1) -> str:
        """Download gbff assemmbly file

        :param acc_id: ncbi accession
        :type acc_id: str
        :param limit: gbff files to return, defaults to 1
        :type limit: int, optional
        :raises ValueError: invalid accession
        :return: file to gbff
        :rtype: str
        """

        search_results, s_count = self.get_search_dict("assembly", acc_id)

        if s_count == 0:
            search_results, s_count = self.get_elink_dict(acc_id)
            if s_count == 0:
                message = f"Query invalid: {acc_id}"
                logger.error(message)
                raise ValueError(message)

        summary_handle = self.entrez.esummary(
            db="assembly",
            rettype="json",
            retmode="json",
            retmax=limit,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )

        # Load summary into dict
        summary = json.loads(summary_handle.read())["result"]

        # Iterate through all results
        gbff_list = []
        for s in summary["uids"][:limit]:
            path = summary[s]["ftppath_genbank"]
            stem = path.split("/")[-1]
            if not stem:
                message = f"Query invalid: {acc_id}"
                logger.error(message)
                raise ValueError(message)
            filepath = f"{stem}_genomic.gbff.gz"
            fna_ftp = f"{path}/{filepath}"
            call(["wget", fna_ftp, "-O", filepath], stderr=subprocess.DEVNULL)
            call(["gunzip", filepath, "-f"])
            gbff_list.append(filepath.rstrip(".gz"))

        return gbff_list[0]

