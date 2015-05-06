"""
Simple script to convert the date-time information given in the particular
Flu alignment into the input date-time file.
"""
import datetime


def str2date_time(instr):
    """
    Convert input string to datetime object.

    Args:
     - instr (str): input string. Accepts one of the formats:
     {MM.DD.YYYY, MM.YYYY, MM/DD/YYYY, MM/YYYY, YYYY}.

    Returns:
     - date (datetime.datetime): parsed date object. If the parsing failed,
     None is returned
    """

    instr = instr.replace('/', '.')
    # import ipdb; ipdb.set_trace()
    try:
        date = datetime.datetime.strptime(instr, "%m.%d.%Y")
    except ValueError:
        date = None
    if date is not None:
        return date

    try:
        date = datetime.datetime.strptime(instr, "%m.%Y")
    except ValueError:
        date = None

    if date is not None:
        return date

    try:
        date = datetime.datetime.strptime(instr, "%Y")
    except ValueError:
        date = None
    return date


def flu_fasta_to_dates():
    """
    Convert fasta file with the flu data into the name,date input csv file.
    Applicable for this given format of the annotation.
    """
    ainf = 'flu.HA.fasta'  # input fasta alignment
    dinf = 'flu.HA.yrs'  # csv dates output

    outstr = []
    aln = AlignIO.read(ainf, 'fasta')
    for a in aln:
        dt = str2date_time(a.name.split('|')[2].strip())
        if dt is not None:
            outstr.append(a.name + ',' +
                     datetime.datetime.strftime(dt, "%Y.%m.%d"))
        with open(dinf, 'w') as outf:
            outf.write('\n'.join(outstr))
    

if __name__ == '__main__':
    pass
