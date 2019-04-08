import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt
import  threedhst.eazyPy   as  eazy
import  threedhst.utils    as  utils


ROOT       = os.environ['SCRATCH']
OUTPUT     = ROOT + '/3D-HST/aegis.v4.1/Eazy'

start      =  1

results    = []
LBGCOLORS  = ['u-g', 'g-r', 'r-i', 'i-z']

for idx in np.arange(41201 - start):
  result = eazy.plotExampleSED(idx= start + idx, pdfname='plot.pdf', MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=OUTPUT,\
                               CACHE_FILE='Same', lrange=[1.e3, 1.e5], axes=None, individual_templates=False, fnu=False)
  
  results.append(result)
            
  if idx % 500 == 0:
    ##  Save as backup during run.
    np.savetxt('dat/aegis_mag_24p5_colors.txt', np.array(results), fmt='%.6lf',\
                                                header='id \t zpeak \t colors \t noisy colors; ' + ''.join(x for x in LBGCOLORS))

##  Save at the run end. 
np.savetxt('dat/aegis_mag_24p5_colors.txt', np.array(results), fmt='%.6lf',\
                                            header='id \t zpeak \t colors \t noisy colors; ' + ''.join(x for x in LBGCOLORS))

print("\n\nDone.\n\n")
