from mdf_connect_client import MDFConnectClient
import mdf_toolbox
import argparse
import time

parser = argparse.ArgumentParser()

#Authentication Args
parser.add_argument("--client_id", type=str, required=True, help="Globus Connect Client ID")
parser.add_argument("--client_secret", type=str, required=True, help="Globus Connect Client Secret")
#Submission Args
parser.add_argument("--title",type=str, help="Dataset title", default="Test MDF Publish 1")  
parser.add_argument("--description", type=str, help="Dataset title", default="Test MDF Publish Description" )  
parser.add_argument("--authors", type=str, help="Dataset Authors", default="ABC")
parser.add_argument("--affiliations", type=str, help="Dataset Author Affiliations", default="UIUC")
parser.add_argument("--related_dois", type=str, help="Related Dataset or Publications", default="xxx.xxxx.xxxxx")
parser.add_argument("--source",type=str, required=True ,help="Data source URL")
parser.add_argument("--update", action="store_true", help="Flag to Update Existing Dataset!", default=False)
parser.add_argument("--test", action="store_true", help="This flags the submission as a test!", default=False)

args = parser.parse_args() 
print(args)

auths = mdf_toolbox.confidential_login(client_id=args.client_id,
                                        client_secret=args.client_secret,
                                        services=["mdf_connect", "mdf_connect_dev"],
                                        make_clients=True)
mdfcc = MDFConnectClient(authorizer=auths['mdf_connect'])
mdfcc.create_dc_block(title=args.title, description=args.description ,authors=args.authors, affiliations=args.affiliations, related_dois=args.related_dois)
mdfcc.add_data_source(args.source)
mdfcc.add_service("mdf_publish")
mdfcc.set_test(args.test)
print("MDF Submission: ", mdfcc.get_submission())
print("MDF Submit Res: ", mdfcc.submit_dataset(update=args.update))
