"""
This module enhances the COSMO-SAC package with the `DirectImport` feature,
providing a more flexible and organized method for managing sigma profiles.
Unlike the standard import method requiring sigma profiles in a single folder,
`DirectImport` allows for storing profiles in separate directories and specifying
the path and name of the sigma file for each component. This approach facilitates
enhanced organization and flexible testing of different sigma profiles for the
same component. The module showcases using `DirectImport` alongside the default
import methods, e.g. DelawareProfileDatabase(), highlighting the new method's
advantages in structuring and accessing sigma profiles.

Author: Ivan Antolovic
E-Mail: ivan.antolovic@tu-berlin.de
"""

import cCOSMO

# Constants
PROFILES_API = "./profiles/pharmaceuticals"
PROFILES_POLY = "./profiles/polymers"
PATH_TO_SIGMAS = "./profiles/sigma3/"
PATH_TO_COMPLIST = "./profiles/complist.txt"
PATH_TO_PROFILES = [PROFILES_API, PROFILES_POLY]

def import_delaware(names, path_to_sigmas, path_to_complist):
    """
    Imports sigma profiles using the DelawareProfileDatabase.

    Parameters:
    - names: A list of component names.
    - path_to_sigmas: The file path to the sigma profiles.
    - path_to_complist: The file path to the component list.

    Returns:
    A COSMO3 object with the imported profiles.
    """
	# Pros/Cons
	# all sigma profiles need to be at the same place
	# all sigma profiles have to be listed in complist.txt
	# best way is to use .xlsx for adding new profiles and convert to .txt
	# most information not really needed to import sigma profile
	# path to sigma profiles + complist.txt
	# names are not equal to identifiers
	# also one name = one identifier --> what if you have modifications?
	# so one name, but different .sigma files
    db = cCOSMO.DelawareProfileDatabase(path_to_complist, path_to_sigmas)
    for name in names:
        iden = db.normalize_identifier(name)
        db.add_profile(iden)
    return cCOSMO.COSMO3(names, db)

def import_direct(names, paths):
    """
    Directly imports sigma profiles using DirectImport.

    Parameters:
    - names: A list of component names.
    - paths: A list of paths to the sigma profiles.

    Returns:
    A COSMO3 object with the imported profiles.
    """
    db = cCOSMO.DirectImport()
    for name, path in zip(names, paths):
        db.add_profile(name, path)
    return cCOSMO.COSMO3(names, db)



if __name__ == "__main__":
    names = ['SIM', 'PLGA50']

    # Example usage of cosmo_delaware and cosmo_direct objects
    cosmo_delaware = import_delaware(names, PATH_TO_SIGMAS, PATH_TO_COMPLIST)
    cosmo_direct = import_direct(names=names, paths=PATH_TO_PROFILES)
