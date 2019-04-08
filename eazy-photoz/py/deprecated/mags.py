>import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt
import  threedhst.eazyPy   as  eazy
import  threedhst.utils    as  utils


ROOT       =  os.environ['SCRATCH']

OUTPUT     =  ROOT + '/MUSE/Rafelski15/EAZY'
## OUTPUT  =  ROOT + '/3D-HST/aegis.v4.1/Eazy'

results    =  []

for idx in np.arange(9970)[::10]:
  ##  eazyPy.py, symlinked in current directory. 
  header, result   = eazy.plotExampleSED(idx=idx, pdfname='plot.pdf', MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY=OUTPUT,\
                                         CACHE_FILE='Same', lrange=[1.e3, 1.e5], axes=None, individual_templates=False, fnu=False)
  
  results.append(result)
            
  if idx % 100 == 0:
    np.savetxt('mags.txt', np.array(results), fmt='%.6lf'.ljust(9), delimiter='\t ', header=header, comments='## ')

np.savetxt('../dat/mags.txt', np.array(results), fmt='%.6lf'.ljust(9), header=header, delimiter='\t ', comments='## ')              

print("\n\nDone.\n\n")
