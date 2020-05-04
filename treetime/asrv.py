import numpy as np
import os
import subprocess
from scipy.stats import gamma
from scipy.special import gammainc


class ASRV(object):
    """
    This class specifies Among-site Rate Variation (ASRV) heterogenous rate calculation
    """
    def __init__(self , alpha=1):
        self.num_catg = 8
        self.alpha = alpha
        self.scale = 1/self.alpha
        
    def parse_alpha(self, logfile):
        """
        parse alpha value for fitted gamma distribution from file
        """
        # path configuration
        path = os.path.abspath(os.path.curdir)
        root_path = path.rpartition('/')[0].replace(" ", "\ ")
        path = os.path.join(root_path, 'data')
        
        proc = subprocess.Popen('grep "^alpha" %s' % (os.path.join(path, logfile)),
                                stdout=subprocess.PIPE, shell=True)
        stdout, stderr = proc.communicate()
        try:
            self.alpha = float(stdout.decode().split(' ')[1])
        except ValueError as val_err:
            print("Failed to parse inferred alpha parameter:", val_err)
            exit()
        except (IndexError, FileNotFoundError) as err:
            print("Failed to open profile file with alpha parameter:", err)
            exit()
            
    
    def calc_rates(self):
        """
        calculate average value in each category as rate in ASRV
        """
        threshold = 1e-6
        
        perc_points = [0]
        perc_points.extend([gamma.ppf(i/self.num_catg, a=self.alpha, scale=self.scale) for i in range(1, self.num_catg)])
        perc_points.append(gamma.ppf(1-threshold, a=self.alpha, scale=self.scale))
        
        rates = np.zeros(self.num_catg)
        
        for i in range(len(perc_points)-1):
            a, b = perc_points[i], perc_points[i+1]
            rates[i] = (gammainc(self.alpha+1, b*self.alpha) - gammainc(self.alpha+1, a*self.alpha)) * self.num_catg
        
        return rates
