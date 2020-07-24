card_options = {}

# Memorly limit
card_options['virtual_memory'] = '100G'
card_options['FORM_compile_cores'] = 10


def process_commands(process_name, process):
    # FORM options
    card_options['FORM_compile_arg'] = process['compile_arg']
    card_options['FORM_compile_optimization'] = process['optimization']
    # QGRAF model
    card_options['qgraf_template_model'] = process['model']
    # Process
    card_options['process_commands'] = ''
    if process['multiparticle'] != '':
        card_options['process_commands'] += 'qgraf_define %(multiparticle)s\n'
    card_options['process_commands'] += 'qgraf_generate %(cmd)s\n'
    card_options['process_commands'] += '!rm -rf TEST_QGRAF_{}_%(order)s\n'.format(
        process_name)
    card_options['process_commands'] += 'output qgraf TEST_QGRAF_{}_%(order)s\n'.format(
        process_name)
    card_options['process_commands'] = card_options['process_commands'] % process


# Processes
processes = {}
processes['epem_a_ddx'] = {
    'cmd': 'e+ e- > a > d d~ /t u s c b QCD^2==0 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttx'] = {
    'cmd': 'e+ e- > a > t t~ /d u s c b QCD^2==0 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ddxg'] = {
    'cmd': 'e+ e- > a > d d~ g/t u s c b QCD^2==1 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxg'] = {
    'cmd': 'e+ e- > a > t t~ g/d u s c b QCD^2==1 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_httx'] = {
    'cmd': 'e+ e- > a > h t t~ /d u s c b QCD^2==0 QED^2==3 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_httxg'] = {
    'cmd': 'e+ e- > a > h t t~ g /d u s c b QCD^2==1 QED^2==3 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxgg'] = {
    'cmd': 'e+ e- > a > t t~ j j/d u s c b QCD^2==2 QED^2==2 []',
    'order': 'LO',
    'multiparticle': 'j = g gh gh~',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxghgh'] = {
    'cmd': 'e+ e- > a > t t~ gh gh~ /d u s c b QCD^2==2 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_httxgg'] = {
    'cmd': 'e+ e- > a > h t t~ j j /d u s c b QCD^2==1 QED^2==3 []',
    'order': 'LO',
    'multiparticle': 'j = g gh gh~',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_httxghgh'] = {
    'cmd': 'e+ e- > a > h t t~ gh gh~ /d u s c b QCD^2==1 QED^2==3 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxggg'] = {
    'cmd': 'e+ e- > a > t t~ j j g/d u s c b QCD^2==3 QED^2==2 []',
    'order': 'LO',
    'multiparticle': 'j = g gh gh~',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxghghg'] = {
    'cmd': 'e+ e- > a > t t~ gh gh~ g/d u s c b QCD^2==3 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxddx'] = {
    'cmd': 'e+ e- > a > t t~ d d~ /u s c b QCD^2==2 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_httxddx'] = {
    'cmd': 'e+ e- > a > h t t~ d d~ /u s c b QCD^2==2 QED^2==3 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}

processes['epem_a_ttxddxg'] = {
    'cmd': 'e+ e- > a > t t~ d d~ g/u s c b QCD^2==3 QED^2==2 []',
    'order': 'LO',
    'multiparticle': '',
    'model': 'epem',
    'optimization': 3,
    'compile_arg': 'all'}


for process_name, process in processes.items():
    with open('qgraf_card.aL', 'r') as card:
        process_commands(process_name, process)
        with open('%s_%s.aL' % (process_name, process['order']), 'w') as f:
            f.write(card.read() % card_options)
