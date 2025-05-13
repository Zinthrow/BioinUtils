from sh import RNAfold
from io import StringIO
import re
from dataclasses import dataclass


@dataclass
class RNAData:
    seq: str
    dot_bracket: str
    dot_count: int
    unpaired_seq_perc: float
    mfe: float


def escape_ansi(line):
    """
    Remove ANSI color coded command line characters from string
    """
    ansi_escape = re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]")
    return ansi_escape.sub("", line)


def format_seq_dict_as_fasta_string_buffer(seq_dict: dict):
    """format seq_dict fasta as string_buffer

    :param seq_dict: dictionary of sequences
    :type seq_dict: dict
    :return: formatted string
    :rtype: str
    """

    stringio = StringIO(
        "\n".join(
            [
                ">" + str(i) + "\n" + seq_dict[i]
                for i in seq_dict
            ]
        )
    )

    return stringio


def input_rfold(fileio, temperature=37):
    """
    Process guides through vienna's RNAfold module
    fileio: filehandle, handle for fasta
    """

    buffer = StringIO()
    RNAfold("--noPS", "-T", temperature, _in=fileio, _out=buffer)

    buffer.seek(0)
    results = escape_ansi(buffer.read()).split(">")[1:]

    r_dict = {}
    for res in results:
        if "\n" not in res:
            print(res)
            continue
        title, seq, dot_notation = res.split("\n")[:-1]
        seq_len = len(seq)
        dot_notation, mfe = dot_notation.split(" (")
        mfe = float(mfe.strip("()"))
        dot_notation = dot_notation[-1 * seq_len:]
        dot_count = dot_notation.count(".")
        dot_perc = round(dot_count / len(dot_notation), 5)

        r_dict[title] = {
            "seq": seq,
            "dot_bracket": dot_notation,
            "dot_count": dot_count,
            "unpaired_seq_perc": dot_perc,
            "mfe": mfe,
        }
    if not r_dict:
        r_dict["NA"] = {
            "seq": None,
            "dot_bracket": None,
            "dot_count": None,
            "unpaired_seq_perc": None,
            "mfe": mfe,
        }

    return r_dict


def convert_dict_to_dataclass(r_dict):
    """
    Convert a dictionary of RNA data into a dataclass.

    :param r_dict: dictionary of results
    :type r_dict: dict
    :return: dict of dataclasses
    :rtype: dict
    """
    dataclass_dict = {}
    
    for title, data in r_dict.items():
        dataclass_instance = RNAData(
            seq=data["seq"],
            dot_bracket=data["dot_bracket"],
            dot_count=data["dot_count"],
            unpaired_seq_perc=data["unpaired_seq_perc"],
            mfe=data["mfe"]
        )
        dataclass_dict[title] = dataclass_instance
    
    return dataclass_dict


def rfold(seq_dict: dict):
    
    seq_io = format_seq_dict_as_fasta_string_buffer(seq_dict)
    r_dict = input_rfold(seq_io)

    return convert_dict_to_dataclass(r_dict)
