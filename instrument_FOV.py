def instrument_FOV(FOV_base, FOV_height):
    
    """

    Modify the file output of Instrument Footprint Editor
    provided by Aladin with a user FOV size.

        http://aladin.u-strasbg.fr/footprint_editor/#
    
    """


    import aladinSAMP
    
    from astropy.io.votable import parse
    
    # download an Aladin Instrument Footprint file from GitHub
    # or making it from http://aladin.u-strasbg.fr/footprint_editor/#
    
    # Reading the Instrument Footprint file template
    votable = parse("footprint_GWsky2.vot")

    # get table in the file (FOV size)
    table = votable.get_first_table()

    
    # get the data in the array member variable:
    data = table.array

    # Convert FOV size degree to arcsecond
    FOV_base_arcsec = FOV_base*3600.0
    FOV_height_arcsec = FOV_height*3600.0

    
    # user FOV size: FOV_base and FOV_height
    data[0] = - FOV_base_arcsec /  2.0,   FOV_height_arcsec / 2.0
    
    data[1] =   FOV_base_arcsec /  2.0,   FOV_height_arcsec / 2.0
    
    data[2] =   FOV_base_arcsec /  2.0, - FOV_height_arcsec / 2.0
    
    data[3] = - FOV_base_arcsec /  2.0, - FOV_height_arcsec / 2.0

    # outputting a VOTable file
    votable.to_xml('instrument_FOV.vot')
    
    # send to aladin the instrument_FOV.vot
    aladinSAMP.send_file('instrument_FOV.vot')

    

     
