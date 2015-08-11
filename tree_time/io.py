import Phylo
import numpy as np
import json
import copy

from tree_time import TreeTime
from tree_anc import TreeAnc, GTR

def _read_json_tree(node, json_clade, data_keys, date_func):
        """
        recursive function to populate tree from the json strcture
        Args:
         - json_clade(dic): the data for the tree node represented as dictionary

         - data_keys(dic): dictionary to convert (some) data keys into the internal 
         tree_time notification

         - date_func(callable): function to convert string date in the json into 
         the datetime object of the tree node
        """
    
        clade_key = 'clade' if 'clade' not in data_keys else  data_keys['clade']
        if clade_key in json_clade:
            node.clade = json_clade[clade_key]

        name_key = 'name' if 'name' not in data_keys else  data_keys['name']
        if name_key in json_clade:
            node.name = json_clade[name_key]

        strain_key = 'strain' if 'strain' not in data_keys else data_keys['strain']
        if strain_key in json_clade:
            node.strain = json_clade[strain_key]

        f_key = 'branch_length' if 'branch_length' not in data_keys else data_keys['branch_length']

        if f_key in json_clade:
            node.branch_length = float(json_clade[f_key])

        
        f_key = 'xvalue' if 'xvalue' not in data_keys else data_keys['xvalue']

        if f_key in json_clade:
            node.xvalue = float(json_clade[f_key])
        
        f_key = 'yvalue' if 'yvalue' not in data_keys else data_keys['yvalue']
        if f_key in json_clade:
            node.yvalue = float(json_clade[f_key])

        f_key = 'days_before_present' if 'days_before_present' not in data_keys else data_keys['days_before_present']
        if f_key in json_clade:
            node.date=float(json_clade[f_key])

        f_key = 'seq' if 'seq' not in data_keys else data_keys['seq']
        if f_key in json_clade:
            node.sequence = np.array(list(json_clade[f_key]))

        f_key = 'yvalue' if 'yvalue' not in data_keys else data_keys['yvalue']
        if f_key in json_clade:
            node.yvalue = float(json_clade[f_key])

        f_key = 'logLH' if 'logLH' not in data_keys else data_keys['logLH']
        if f_key in json_clade:
            node.logLH = float(json_clade[f_key])

        if len(node.clades):
            json["children"] = []
            for ch in node.clades:
                json["children"].append(self.to_json(ch))

def from_json(cls, inf, json_keys={"branch_len":"xvalue"}, date_func=lambda x: None):
    """
    Load tree from json file. 

    Args:
     - inf(str): pth to the input file
     - json_keys(dic): names of the parameters in the json file. The names 
     for the following parameters should be specified: 
     {"date": ..., "branch_len": ...}
     - date_func(callable): function to convert the date representation from 
     json parameter string to datetime object (will be assigned as raw_date)
    Returns:
     - TreeTime object 
    """
    with open (inf) as json_file:
        data = json.load(json_file)
    
    if len(data) < 1 or 'children' not in data:
        raise IOError("Wrong format of json file")
    t = Phylo.BaseTree.Tree()
    ttime = TreeTime(t)

    ttime.read_json_tree(t.root, data, json_keys, date_func)
           
    return ttime


def _read_dates_file(self, inf, **kwargs):
        """
        Read dates from the file into python dictionary. The input file should
        be in csv format 'node name, date'. The date will be converted to the
        datetime object and added to the dictionary {node name: datetime}

        Args:
         - inf(str): path to input file

        KWargs:
         - verbose(int): how verbose should the output be

        Returns:
         - dic(dic): dictionary  {NodeName: Date as datetime object}
        """

        def str_to_date(instr):
            """
            Convert input string to datetime object.
    
            Args:
             - instr (str): input string. Accepts one of the formats:
             {YYYY.MM.DD, YYYY.MM, YYYY}.
    
            Returns:
             - date (datetime.datetime): parsed date object. If the parsing
             failed, None is returned
            """
            try:
                date = datetime.datetime.strptime(instr, "%Y.%m.%d")
            except ValueError:
                date = None
            if date is not None:
                return date
    
            try:
                date = datetime.datetime.strptime(instr, "%Y.%m")
            except ValueError:
                date = None
    
            if date is not None:
                return date
    
            try:
                date = datetime.datetime.strptime(instr, "%Y")
            except ValueError:
                date = None
    
            return date

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = 10

        if verbose > 3:
            print ("Reading datetime info for the tree nodes...")
        with open(inf, 'r') as finf:
            all_ss = finf.readlines()
        if verbose > 5:
            print ("Loaded %d lines form dates file" % len(all_ss))
        try:
            dic = {s.split(',')[0]: str_to_date(s.split(',')[1].strip())
                   for s in all_ss if not s.startswith("#")}
            if verbose > 3:
                print ("Parsed data in %d lines of %d input, %d corrupted"
                       % (len(dic), len(all_ss), len(all_ss) - len(dic)))
            return dic
        except ValueError:
            # unable to read all dates, the file is corrupted - go one by one
            print ("Unable to perform parsing of the dates file, file is "
                   "corrupted. Return empty dictionary.")
            return {}

if __name__=='__main__':
    pass
