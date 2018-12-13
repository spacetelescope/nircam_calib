"""

Create YAML configuration file to accompany the flight2ground.py script.

Authors:
--------
    - Alicia Canipe

Dependencies:
-------------

    - yaml

Use:
----
    python create_config.py --output_name='flight2ground_config'

"""

# Required packages
import argparse
import yaml


def main(args):
    """Main function.

    """

    # Create the dictionary with mnemonics and headers to add
    config = dict(mnemonics=dict(IFGS_GS1='IFGS_GUIDERSTATE1',
                                 IFGS_GS2='IFGS_GUIDERSTATE2',
                                 IFGS_ES='IFGS_EXPOSE_STAT',
                                 IMIR_ES='IMIR_EXP_STAT',
                                 IMIR_FS='IMIR_FILTER_STAT',
                                 IMIR_GS='IMIR_GRATING_STAT',
                                 IMIRFWCP='IMIR_HK_FW_CUR_POS'
                                 ),
                  headers=dict(COLCORNR=1,
                               ROWCORNR=1,
                               COLSTART=1,
                               COLSTOP=1,
                               )
                  )

    # Create the YAML file and dump the dictionary into it
    if args.output_name is not 'None':
        with open(args.output_name+'.yml', 'w') as configfile:
            yaml.dump(config, configfile, default_flow_style=False)
    else:
        with open('flight2ground_configs.yml', 'w') as configfile:
            yaml.dump(config, configfile, default_flow_style=False)


if __name__ == '__main__':

    # Command line argument handler.
    parser = argparse.ArgumentParser(
        description='Create YAML configuration file to accompany the flight2ground.py script.',
        epilog='example: create_config.py --output_name')
    parser.add_argument('--output_name', help='config file name', default='None')

    args = parser.parse_args()
    main(args)
