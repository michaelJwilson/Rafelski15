##  import  matplotlib;                matplotlib.use('AGG')

import  copy
import  json
import  glob
import  numpy              as      np
import  pylab              as      pl

import  matplotlib.cm      as      cm
import  matplotlib.pyplot  as      plt
import  pandas             as      pd

from    collections        import  OrderedDict
from    pylab              import  rcParams
from    app_mags           import  get_colors
from    utils              import  latexify
from    colourcut_dNdz     import  colourcut
from    astropy.table      import  Table
from    depths             import  get_depths


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

degraded_depths = get_depths()

np.random.seed(seed=314)

def colourtrack(ttype = 'g'):    
    '''
    Plot locus of colour box for given dropout selection. 
    '''
                  ##  Hildebrandt (2018), https://arxiv.org/pdf/0903.3951.pdf;  Otherwise, GoldRush (Ono, 2018).  
                  ##  u-drops:  1.5 < (u-g) and -1 < (g-r) < 1.2 and 1.5 (g-r) < (u-g) - 0.75.
                  ##      BzK:  

    selections = {'u':   {'bcol': 'u-g', 'rcol': 'g-r', 'minbcol': 1.5, 'maxrcol': 1.2, 'gradient': 1.5, 'intercept': 0.75, 'hiz': 4.5},\
                  'g':   {'bcol': 'g-r', 'rcol': 'r-i', 'minbcol': 1.0, 'maxrcol': 1.0, 'gradient': 1.0, 'intercept': 1.00, 'hiz': 6.0},\
                  'r':   {'bcol': 'r-i', 'rcol': 'i-z', 'minbcol': 1.2, 'maxrcol': 0.7, 'gradient': 1.5, 'intercept': 1.00, 'hiz': 7.5},\
                  'z':   {'bcol': 'i-z', 'rcol': 'z-y', 'minbcol': 1.5, 'maxrcol': 0.5, 'gradient': 2.0, 'intercept': 1.10, 'hiz': 7.5}}
            
    if ttype == 'BzK':
        bcols      =  np.arange(-4.0, 7.0, 0.001)
        rcols      =  bcols -0.2

        pl.plot(rcols, bcols, 'k-')
        pl.plot(2.5 * np.ones_like(rcols[rcols > 2.5]), bcols[rcols > 2.5], 'k--')

    elif ttype == 'Euclid':
        pass

    elif ttype == 'Synergy':
        pass

    else:
        ##  Dropouts.
        bcol       =  selections[ttype]['bcol']
        rcol       =  selections[ttype]['rcol']
        
        minbcol    =  selections[ttype]['minbcol']        ## Detect the break
        maxrcol    =  selections[ttype]['maxrcol']        ## Flat spectra above break. 
        
        gradient   =  selections[ttype]['gradient']
        intercept  =  selections[ttype]['intercept']

        rcols      =  np.arange(-5.0, 10.0, 0.001)
        bcols      =  gradient * rcols + intercept

        minrcol    = (minbcol - intercept) / gradient
        bcollim    =  gradient * maxrcol + intercept

        pl.plot(rcols[rcols < minrcol], minbcol * np.ones_like(rcols[rcols < minrcol]), 'k', lw=0.6)
        pl.plot(maxrcol * np.ones_like(bcols[bcols > bcollim]), bcols[bcols > bcollim], 'k', lw=0.6)
    
        ##  Gradient
        pl.plot(rcols[(rcols > minrcol) & (bcols < bcollim)], bcols[(rcols > minrcol) & (bcols < bcollim)], c='k', linestyle='-', lw=0.6)

        if ttype == 'u':
            ##  Hildebrandt u-drops has a lower limit on (g-r)
            pl.plot(-1. * np.ones_like(bcols[bcols > minbcol]), bcols[bcols > minbcol], c='k', linestyle='--', lw=0.6)
            
def plot_cat(ttype = 'g', DODEPTH='FULL'):    
    if   ttype == 'u':
         bcol = 'u-g'
         rcol = 'g-r' 

         bmin = -0.5
         bmax =  2.5 

         rmin = -0.3 
         rmax =  1.2 

         infrac  = 1.00
         outfrac = 0.10
         
    elif ttype == 'g':
         bcol = 'g-r'
         rcol = 'r-i'
        
         bmin = -0.5 
         bmax =  2.5 

         rmin = -0.3 
         rmax =  1.2 

         infrac  = 1.00
         outfrac = 0.10

    elif ttype == 'BMBX':
         bcol = 'g-r'
         rcol = 'r-i'

         bmin = -0.5 
         bmax =  2.5

         rmin = -0.3 
         rmax =  1.2 

         infrac  = 1.00
         outfrac = 0.10

    elif ttype == 'BzK':
         bcol = 'B-z'
         rcol = 'z-K'

         bmin = -1.5 
         bmax =  4.0 

         rmin = -1.0 
         rmax =  7.0 

         infrac  = 0.20
         outfrac = 0.10

    elif ttype == 'Euclid':
         bcol = 'Y-J'
         rcol = 'J-H'

         if DODEPTH == 'Full':
             bmin = -0.10
             bmax =  1.00

             rmin = -0.25
             rmax =  1.60

             infrac  = 0.50
             outfrac = 0.10

         else:
             bmin = -2.5
             bmax =  1.5

             rmin = -2.5
             rmax =  1.5

             infrac  = 1.00
             outfrac = 0.10

    elif ttype == 'Synergy':
         bcol = 'r-i'
         rcol = 'i-H'

         bmin = -1.0
         bmax =  4.0

         rmin = -0.5
         rmax =  6.0

         infrac  = 1.00
         outfrac = 0.10

    else:
        raise ValueError("\n\nSelection type %s is not available with cat_tracks.\n\n" % ttype)

    ##  -- Start -- 
    files         =  glob.glob('colors/scratch/*.fits')
    files         = [x.split('.fits')[0] for x in files]

    files         = [x.split('colors/scratch/3DHST_')[1] for x in files]

    fields        = [x.split('_')[0] for x in files]
    ufields       = list(set(fields))
    nfields       = len(ufields)

    depths        = [xx.split('_')[1].split('_')[0] for xx in files]
    depths        = ['Full' if depth != 'Degraded' else depth for depth in depths]
    udepths       =  set(depths)

    count         =  0

    detected      =  0
    udetected     =  0

    for ii, xx in enumerate(files):
        depth  = depths[ii]
        field  = fields[ii]

        drops  = []

        if (field == 'UVUDF') & (depth == DODEPTH):
            fname = 'colors/scratch/3DHST_' + xx + '.fits'

            print('\n\nSolving for %s' % fname)

            magt  = Table.read(fname, format='fits')
            cols  = magt.colnames

            mags  = copy.copy(cols)

            mags.remove('id')
            mags.remove('zpeak')

            count += 1

            for row in magt:
              if (row['zpeak'] > -99.):
                rowmags    = [row[x] for x in mags]
                magdict    = dict(zip(mags, rowmags))

                ##  'z-K'
                colors     = get_colors(magdict, get_colors=['g-r', 'r-i', 'u-g', 'g-r', 'u-z', 'B-z', 'J-H', 'Y-J', 'g-i', 'i-H'], fname = None)

                if ttype     == 'u':
                  is_detected = (magdict['r'] < degraded_depths['r'])

                elif ttype == 'g':
                  is_detected = (magdict['i'] < degraded_depths['i'])

                elif ttype == 'BzK':
                  is_detected = (magdict['K'] < degraded_depths['K'])  

                elif ttype in ['Euclid', 'Synergy']:
                  is_detected = (magdict['H'] < degraded_depths['H'])  

                else:
                  is_detected = False

                if DODEPTH == 'Full':
                  is_detected =  True

                if is_detected:
                  detected   += 1  
                  is_lbg      = colourcut(magdict, dropband=ttype, good=True, fourthlimit=False, BzK_type='all')

                  if ttype == 'BzK':
                      BzK     = colors['z-K'] - colors['B-z']

                      print('%+.4lf \t %+.4lf \t %+.4lf \t %+.4lf \t %s' % (magdict['B'], magdict['z'], magdict['K'], BzK, str(is_lbg)))

                  else:
                      print('%+.4lf \t %+.4lf \t %+.3lf \t %+.3lf \t %s \t %+.3lf \t %+.3lf' % (magdict['u'], magdict['g'], magdict['r'],\
                                                                                                row['zpeak'], str(is_lbg), degraded_depths['r'],\
                                                                                                degraded_depths['i']))

                  draw = np.random.uniform(0.0, 1.0, 1)

                  if ttype == 'Euclid':
                    if   is_lbg[0] & (draw <= infrac):
                      cax  = plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], s=10, vmin=0.0, vmax=5.0, rasterized=True, alpha=1.)

                    elif is_lbg[1] & (draw <= infrac):
                      cax  = plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], s=20, vmin=0.0, vmax=5.0, rasterized=True, alpha=1., marker='*')  

                    else:
                      if (draw <= outfrac):
                          plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], marker='x', alpha=0.2, s=10, vmin=0.0, vmax=5.0, rasterized=True)

                  else:
                    if is_lbg & (draw <= infrac):
                      cax  = plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], s=10, vmin=0.0, vmax=5.0, rasterized=True, alpha=1.)

                    else:
                      if (draw <= outfrac):                          
                          plt.scatter(colors[rcol], colors[bcol], c=row['zpeak'], marker='x', alpha=0.2, s=10, vmin=0.0, vmax=5.0, rasterized=True)

                else:
                    udetected += 1

                    ##  print('NO DETECTION:  %+.4lf \t %+.4lf \t %+.3lf \t %+.3lf \t %+.3lf \t %+.3lf' % (magdict['u'], magdict['g'], magdict['r'],\
                    ##                                                                                     row['zpeak'], degraded_depths['r'],\
                    ##                                                                                     degraded_depths['i']))

            print('Detection rates: %d \t %d' % (detected, udetected))

    ax  = pl.gca()
    fig = pl.gcf()

    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')

    pl.xlabel(r'$%s$' % rcol)
    pl.ylabel(r'$%s$' % bcol)
    
    stretch = 0.0

    pl.xlim(rmin - stretch, rmax + stretch)
    pl.ylim(bmin - stretch, bmax + stretch)

    ##  [left, bottom, width, height].
    cbaxes = fig.add_axes()   ## [0.8, 0.1, 0.03, 0.8]
    cb     = plt.colorbar(ax=ax, cax=cbaxes, label=r'redshift')  

    cb.set_alpha(1.0)
    cb.draw_all()

    if DODEPTH == 'Full':
        cb.set_alpha(0.0)
        cb.set_label('')
        cb.set_ticklabels([])
        cb.set_ticks([])
        cb.draw_all()

    return  bcol, rcol, [[rmin, rmax], [bmin, bmax]]


if __name__ == "__main__":
  print("\n\nWelcome to colour tracks.\n\n")

  for DODEPTH in ['Full']:      ##  ['Full', 'Degraded']
    for ttype in ['Euclid']:    ##  ['u', 'g', 'BzK', 'Euclid', 'Synergy']
      pl.clf()
      
      bcol, rcol, [[rmin, rmax], [bmin, bmax]] = plot_cat(ttype=ttype, DODEPTH = DODEPTH)  

      colourtrack(ttype=ttype)

      plt.tight_layout()     

      pl.show()
      ##  pl.savefig('plots/ccplots/%s_colour_track_%s.png' % (ttype, DODEPTH))
      
  print("\n\nDone.\n\n")
