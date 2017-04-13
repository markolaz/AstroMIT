def testtypping(cfgfile="example.cfg"):
    infile = "LC/114913.tesslc"
    lc = Lightcurve().create_lightcurve(infile,cfgfile=cfgfile)
    lc.gen_jd()
    lc.plot_lc()
    infile = "LC/kplr010000069-2009166043257_llc.fits"
    lc = Lightcurve().create_lightcurve(infile,cfgfile=cfgfile)
    lc.load_from_file()
    lc.plot_lc()
    detrend_func = COSDetrend(cfgfile=cfgfile)
    lc.detrend_lc(detrend_func)
    lc.plot_lc(label='ltflc')
    lc.write_to_file("LC/kplr010000069-2009166043257_ltf.lc")

    return

def test(cfgfile="example.cfg"):
    infile = "LC/114913.tesslc"
    lc = AsciiLightcurve(cols={'jd': 2, 'rlc': 8, 'ltflc': 12,
                               'x': 3, 'y': 4, "bg": 12, "bgerr": 12, 'cadence': 2})
    lc.name = infile
    lc.gen_jd()
    lc.plot_lc()
    detrend_func = FocusDetrend(cfgfile=cfgfile)
    lc.detrend_lc(detrend_func)
    lc.plot_lc(label='ltflc')
    detrend_func = COSDetrend(cfgfile=cfgfile)
    lc.detrend_lc(detrend_func)
    lc.plot_lc(label='ltflc')
    detrend_func = PolyDetrend(cfgfile=cfgfile)
    lc.detrend_lc(detrend_func)
    lc.plot_lc(label='ltflc')
    lc.write_to_file("LC/114913.ltf")
    return


