import os

def write_ispyb_xml(filename, full_command_line, write_directory, xml_results):
    '''Write items in the _xml_results into an XML file to be stored in
    ISPyB'''
    xml_template = os.path.join(os.environ['FAST_EP_ROOT'],
                                'lib', 'templates', 'ispyb.xml')
    phs_stat_fom_template = os.path.join(
        os.environ['FAST_EP_ROOT'], 'lib', 'templates',
        'phasing_statistics_fom.xml')
    phs_stat_mapcc_template = os.path.join(
        os.environ['FAST_EP_ROOT'], 'lib', 'templates',
        'phasing_statistics_mapcc.xml')

    if not os.path.exists(xml_template):
        print 'XML template not found: %s' % xml_template
        return
    if not os.path.exists(phs_stat_fom_template):
        print 'XML template not found: %s' % phs_stat_fom_template
        return
    if not os.path.exists(phs_stat_mapcc_template):
        print 'XML template not found: %s' % phs_stat_mapcc_template
        return

    # get phasing statistics from xml_results

    (all_phs_stat_fom, all_phs_stat_mapcc) = get_phasing_statistics(
        phs_stat_fom_template, phs_stat_mapcc_template, xml_results)

    import datetime
    time_stamp = '%4d-%02d-%02d %02d:%02d:%02d' % tuple(datetime.datetime.now(
        ).timetuple()[:6])
    
    open(filename, 'w').write(
        open(xml_template, 'r').read().format(
            commandline = full_command_line,
            results_directory = write_directory,
            spacegroup_id = xml_results['SPACEGROUP'],
            solvent_content = xml_results['SOLVENTCONTENT'],
            enantiomorph = xml_results['ENANTIOMORPH'],
            lowres = xml_results['LOWRES'],
            highres = xml_results['HIGHRES'],
            shelxc_spacegroup = xml_results['SHELXC_SPACEGROUP_ID'],
            phasing_statistics_fom = all_phs_stat_fom,
            phasing_statistics_mapcc = all_phs_stat_mapcc,
            time_stamp = time_stamp
            ))

def get_phasing_statistics(fom_template, cc_template, xml_results):
    total_bins = 1
    all_phs_stat_fom = ""
    all_phs_stat_mapcc = ""

    # find number of bins - use RESOLUTION_LOW as the field to check for this

    done = False
    while not done:
        bin_number_name = str(total_bins).zfill(2)
        try:
            resolution_low = xml_results['RESOLUTION_LOW' + bin_number_name]
        except KeyError:
            done = True
            continue
        total_bins += 1
        
    for bin_number in xrange(total_bins):
        bin_number_name = str(bin_number).zfill(2)
        resolution_low = float(xml_results['RESOLUTION_LOW' + bin_number_name])
        resolution_high = float(xml_results['RESOLUTION_HIGH' + bin_number_name])
        fom = float(xml_results['FOM' + bin_number_name])
        mapcc = float(xml_results['MAPCC' + bin_number_name])
        nreflections = int(xml_results['NREFLECTIONS' + bin_number_name])
        if resolution_low == None or resolution_high == None or \
          fom == None or mapcc == None or nreflections == None:
            raise RuntimeError, "One of the fields is empty."
        all_phs_stat_fom += open(fom_template,'r').read().format(
                bin_number = bin_number + 1, 
                number_bins = total_bins,
                bin_low_res = resolution_low,
                bin_high_res = resolution_high,
                bin_fom = fom,
                num_refl = nreflections)
        all_phs_stat_mapcc +=  open(cc_template,'r').read().format(
                bin_number = bin_number + 1,
                number_bins = total_bins,
                bin_low_res = resolution_low,
                bin_high_res = resolution_high,
                bin_map_cc = mapcc,
                num_refl = nreflections)
    return (all_phs_stat_fom, all_phs_stat_mapcc)
        
