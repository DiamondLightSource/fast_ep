import os

def write_ispyb_xml(filename, xml_results):
    '''Write items in the _xml_results into an XML file to be stored in ISPyB'''
    xml_template = os.path.join(os.environ['FAST_DP_ROOT'],
                                'lib', 'templates', 'ispyb.xml')
    phasing_statistics_fom_template = os.path.join(os.environ['FAST_DP_ROOT'],
                                'lib', 'templates', 'phasing_statistics_fom.xml')
    phasing_statistics_mapcc_template = os.path.join(os.environ['FAST_DP_ROOT'],
                                'lib', 'templates', 'phasing_statistics_mapcc.xml')

    if not os.path.exists(xml_template):
        print 'XML template not found: %s' % xml_template
        return
    if not os.path.exists(phasing_statistics_fom_template):
        print 'XML template not found: %s' % phasing_statistics_fom_template
        return
    if not os.path.exists(phasing_statistics_mapcc_template):
        print 'XML template not found: %s' % phasing_statistics_mapcc_template
        return

    open(filename, 'w').write(
        open(xml_template, 'r').read().format(
            spacegroup = xml_results['SPACEGROUP']
            ))