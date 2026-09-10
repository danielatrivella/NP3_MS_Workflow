## INSTALL
# use conda or mamba to install the environment
conda env create -f environment_np3_py4cy.yml 

# activate the environment
conda activate np3_py4cy

## USE
# check the network mounting script arguments - help usage
python network_mount_py4cy.py -h
# execute the network mounting script with cytoscape opened in the desired or an empty project
# the argument -y is optional
python network_mount_py4cy.py -s <ssmn_selfloops> -p <ssmn_protonated_selfloops> -c <clean_table_file_csv> -y <cytoscape_styles_xml>

