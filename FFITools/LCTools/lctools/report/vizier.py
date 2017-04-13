import numpy as np
import pyfits
import urllib
from astropy.io import votable
from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u

def tic_query(ticid):
    """
    query the star's information from tic database on tesslate, 
    result is an array of dict, each dict is a row from the database. For row content, see tsig wiki TIC3 Schema,
    ticid is required to be a string.
    """
    import tsig
    from tsig import catalog
    c = tsig.catalog.TIC()
    result = c.query_by_id([ticid])
    TICheader = open("info/ticheader.txt").readline() 
    dic = {}
    for key,value in zip(TICheader.split(','), result[0]):
        dic[key.strip()] = value
    return dic

def sed_from_tic(starid, ra, dec):
    queryresult = tic_query(starid)
    # print queryresult
    filterinfo = np.recfromcsv("info/standardfilters.csv", delimiter=',')
    filter_dict = {'KMAG': '2MASS Ks', 'JMAG': '2MASS J', 'HMAG': '2MASS H', 'W1MAG': 'WISE-1',
                   'W2MAG': 'WISE-2', 'UMAG': 'SDSS u', 'GMAG': 'SDSS g', \
                           'RMAG': 'SDSS r', 'IMAG': 'SDSS i', 'ZMAG': 'SDSS z', 'BMAG': 'Johnson:B', 'VMAG': 'Johnson:V'}
    filter_all = filterinfo['band']
    wavelength_all = filterinfo['lambda']
    dlambda_all = filterinfo['dlambda']
    Z0_all = filterinfo['zeropoint']
    mags=[]; e_mags=[]; filters=[]; wavelength=[]; dlambda=[]; Z0=[]
    for key in filter_dict.keys():
        if queryresult[key] == 0.0:
            continue
        ind = key == filter_all
        filters.append(key)
        wavelength.append(wavelength_all[ind][0])
        dlambda.append(dlambda_all[ind][0])
        Z0.append(Z0_all[ind][0])
        mags.append(queryresult[key])
        e_mags.append(queryresult['E_'+key])
    
    mags=np.array(mags)
    e_mags=np.array(e_mags)
    wavelength = np.array(wavelength)
    dlambda = np.array(dlambda)
    Z0 = np.array(Z0)
    #print len(mags), len(Z0)
    #flux is in erg/s/cm^2/Angstrom 
    flux = 10**(-(mags+Z0)/2.5)
    eflux = flux*e_mags 
    #use band width and band centre to convert to Jansky

    fluxmv = (flux/3.e-5)*(wavelength**2.-dlambda**2.)
    efluxmv = (eflux/3.e-5)*(wavelength**2.-dlambda**2.)
    #print filters
    #print wavelength
    c = coord.SkyCoord(float(ra) * u.deg, float(dec) * u.deg, frame='icrs')
    tb = IrsaDust.get_extinction_table(c)
    filters2 = tb['Filter_name']
    allA = tb['A_SandF']
    A = []
    f_ob_orig = []
    wav_ob = []
    unc = []
    
    for i in xrange(len(filters)):
        
        ind = np.where(filters2 == filter_dict[filters[i]])[0]
        if len(ind) == 0:
            continue
        A.append(np.mean(allA[ind]))
        
        f_ob_orig.append(np.mean(fluxmv[i]))
        wav_ob.append(np.mean(wavelength[i]))
        unc.append(np.mean(efluxmv[i]))
    
    return [A, f_ob_orig, wav_ob, unc]

def sed_from_web(starid, ra, dec):
    urllib.urlretrieve('http://vizier.u-strasbg.fr/viz-bin/sed?-c=' + ra + "%2C" + dec + '&-c.rs=1', starid + '_sed.vot')
    print 'http://vizier.u-strasbg.fr/viz-bin/sed?-c=' + ra + "%2C" + dec + '&-c.rs=0.005'

    tb = votable.parse_single_table(starid + '_sed.vot')
    data = tb.array
    wav_all = 3e5 * 1e4 / data['_sed_freq'].data  # angstrom
    f_all = data['_sed_flux'].data
    unc_all = data['_sed_eflux'].data
    filters = data['_sed_filter'].data

    #print filters
    filter_dict = {'2MASS:Ks': '2MASS Ks', '2MASS:J': '2MASS J', '2MASS:H': '2MASS H', 'WISE:W1': 'WISE-1',
                   'WISE:W2': 'WISE-2', 'SDSS:u': 'SDSS u', 'SDSS:g': 'SDSS g', \
                   'SDSS:r': 'SDSS r', 'SDSS:i': 'SDSS i', 'SDSS:z': 'SDSS z'}

    c = coord.SkyCoord(float(ra) * u.deg, float(dec) * u.deg, frame='icrs')
    tb = IrsaDust.get_extinction_table(c)
    filters2 = tb['Filter_name']
    allA = tb['A_SandF']
    A = []
    f_ob_orig = []
    wav_ob = []
    unc = []

    for f in filters:
        if f in filter_dict.keys():
            filtmatch = filter_dict[f]
            ind = np.where(filters2 == filtmatch)[0]
            A.append(np.mean(allA[ind]))
            ind = np.where(filters == f)[0]
            f_ob_orig.append(np.mean(f_all[ind]))
            wav_ob.append(np.mean(wav_all[ind]))
            unc.append(np.mean(unc_all[ind]))
    return [A,f_ob_orig, wav_ob, unc]

def test_sed(starid, ra, dec):

    A1, f_ob_orig1, wav_ob1, unc1 = sed_from_web(starid, ra, dec)
    # print np.unique(np.array(A1))
    # print np.unique(np.array(f_ob_orig1))
    # print np.unique(np.array(wav_ob1))
    # print np.unique(np.array(unc1))
    
    A2, f_ob_orig2, wav_ob2, unc2 =  sed_from_tic(starid, ra, dec) 
    # print A2
    # print f_ob_orig2
    # print wav_ob2
    # print unc2
    return

def sed(starid, ra, dec, ax):
    """
    Retrieves star's magnitudes in different bands from ExoFOP and fits SED from Castelli & Kurucz library.
    Note: you need an account to access Kepler ExoFOP data. I do not have one, so I'm using the format for K2 targets
    instead.

    :param kic: K2 EPIC ID of target.
    :param ra: RA of target.
    :param dec: Dec of target.
    :param ax: plot handle of SED plot.

    :return: ax: SED plot.
    """
    try:
        A,f_ob_orig, wav_ob, unc =  sed_from_tic(starid, ra, dec) 
    except:
        A,f_ob_orig, wav_ob, unc = sed_from_web(starid, ra, dec)
    f_ob_orig = np.array(f_ob_orig)
    A = np.array(A)
    wav_ob = np.array(wav_ob)
    f_ob = (f_ob_orig * 10 ** (A / 2.5))

    metallicity = ['ckp00']
    m = [0.0]
    t = np.arange(3500, 13000, 250)
    t2 = np.arange(14000, 50000, 1000)
    t = np.concatenate((t, t2))

    log_g = ['g20', 'g25', 'g30', 'g35', 'g40', 'g45', 'g50']
    g = np.arange(2., 5., 0.5)

    best_m = 0
    best_t = 0
    best_g = 0
    best_off = 0.0
    chi2_rec = 1e6

    # Do grid search to find best-fit SED
    # This loop can probably be parallelized to save time
    for im, mval in enumerate(m):
        for it, tval in enumerate(t):
            for ig, gval in enumerate(g):
                # load model
                hdulist = pyfits.open('SEDdata/' + metallicity[im] + '/' + metallicity[im] + '_' + str(tval) + '.fits')
                data = hdulist[1].data
                wmod = data['WAVELENGTH']
                fmod = data[log_g[ig]] * 3.34e4 * wmod ** 2

                # fit observations
                f_int = np.exp(np.interp(np.log(wav_ob), np.log(wmod), np.log(fmod)))
                offsets = np.linspace(np.log(min(f_ob / f_int)), np.log(max(f_ob / f_int)), 51)
                for i_off, offset in enumerate(offsets):
                    chi2 = sum((f_int * np.exp(offset) - f_ob) ** 2)

                    # print 'chi2=', chi2, mval, tval, gval
                    if chi2 < chi2_rec:
                        chi2_rec = chi2
                        best_m = im
                        best_g = ig
                        best_t = it
                        best_off = offset

    print 'best fit: m=', m[best_m], 'T=', t[best_t], 'log g=', g[best_g]

    hdulist = pyfits.open('SEDdata/' + metallicity[best_m] + '/' + metallicity[best_m] + '_' + str(t[best_t]) + '.fits')
    data = hdulist[1].data
    wmod = data['WAVELENGTH']
    fmod = data[log_g[best_g]] * 3.34e4 * wmod ** 2
    fmod *= np.exp(best_off)

    ax.plot(wmod / 1e4, fmod, label='Castelli & Kurucz model')
    ax.set_xscale('log')
    ax.plot(wav_ob / 1e4, f_ob_orig, lw=0, marker='s', label='Uncorrected', ms=10)
    ax.plot(wav_ob / 1e4, f_ob, lw=0, marker='o', label='Corrected for extinction', ms=10)
    ax.set_xlabel(r'${\rm Wavelength} \ (\mu m)}$', fontsize=18)
    ax.set_xlim(0.1, max(wmod) / 1e4)
    ax.set_ylabel(r'$F_{\nu} \ {\rm (Jy)}$', fontsize=18)
    ax.legend()
    return ax
