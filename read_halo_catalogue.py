""" Modifies the dragons.nbody.read_halo_catalog function to not read in all
halos at a given snapshot.
"""
from os import listdir
from os import path

#Set data type for all of the halo properties in the catalog
catalog_halo_dtype = np.dtype(dict(names=("id_MBP", "M_vir", "n_particles",
                                          "position_COM", "position_MBP",
                                          "velocity_COM", "velocity_MBP",
                                          "R_vir", "R_halo", "R_max", "V_max",
                                          "sigma_v", "spin", "q_triaxial",
                                          "s_triaxial", "shape_eigen_vectors",
                                          "padding"),
                                   formats=('q', 'f8', 'i4', ('f4', 3),
                                            ('f4', 3), ('f4', 3), ('f4', 3),
                                            'f4', 'f4', 'f4', 'f4', 'f4',
                                            ('f4', 3), 'f4', 'f4',
                                            ('f4', (3, 3)), 'S8')),
                              align=True)

#Set data type for header properties in the catalog
catalog_header_dtype = np.dtype(dict(names=("i_file", "N_files",
                                            "N_halos_file", "N_halos_total"),
                                     formats=['i4', ]*4), align=True)

#Function to read the catalog at location catalog_loc
def read_halo_catalog(catalog_loc):

    """ Read in a halo catalog produced by gbpCode.

    *Args*:
        catalog_loc : str
            Full path to input catalog file or directory

    *Returns*:
        halo : array
            The catalog of halos
    """

    #Converts str to list
    if type(catalog_loc) is str:
        catalog_loc = [catalog_loc, ]

    #If catalog_loc is an existing directory
    if path.isdir(catalog_loc[0]):
        dirname = catalog_loc[0] #Directory name (str)
        catalog_loc = listdir(catalog_loc[0]) #List of all files in the directory
        file_num = [int(f.rsplit(".")[-1]) for f in catalog_loc] #lists numbers
        #at the end of files in the directory, e.g. 48 for
        #subfind_019.catalog_groups_properties.48. Seem to range from 0 to 63
        #Sort the files by their file_num:
        sort_index = np.argsort(file_num)
        catalog_loc = [catalog_loc[i] for i in sort_index]
        #Create list of paths to each file in the directory/catalogue
        catalog_loc = [path.join(dirname, f) for f in catalog_loc]
    n_halos = np.fromfile(catalog_loc[0], catalog_header_dtype,
                          1)[0]["N_halos_total"] #Construct an array containing
    #the header from the first file, select the 0th array which contains the data
    #values, and select the value corresponding to "N_halos_total"
    halo = np.empty(n_halos, dtype=catalog_halo_dtype) #create an empty array of
    #shape (n_halos,)
    print "Reading in {:d} halos...".format(n_halos)

    n_halos = 0
    for f in tqdm(catalog_loc): #Loops through the catalogue, tqdm makes a
    #progress bar
        with open(f, "rb") as fd: #open each file
            n_halos_file = \
                np.fromfile(fd, catalog_header_dtype, 1)[0]["N_halos_file"]
                #Load the number of halos within each file
            halo[n_halos:n_halos+n_halos_file] = \
                np.fromfile(fd, catalog_halo_dtype, n_halos_file) # Load in the
                #n_halos_file elements in the file
        n_halos += n_halos_file
                #increment the number of halos count by n_halos_file

    return halo[list(catalog_halo_dtype.names[:-1])]
