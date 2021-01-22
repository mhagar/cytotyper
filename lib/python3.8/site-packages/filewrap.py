#!/usr/bin/env python3

import os, tarfile, gzip, zipfile, zlib

# Copy target directory and all of it's subdirectories/files to a destination directory.
def copydir(destination_path, target_path):
    if (type(destination_path) is not str):
        raise TypeError("The given file path " + str(destination_path) + " isn't a string!")

    if (type(target_path) is not str):
        raise TypeError("The given file path " + str(target_path) + " isn't a string!")

    if (path_exists(destination_path) == False):
        raise FileNotFoundError("The given file path " + str(destination_path) + " doesn't exist!")

    if (path_exists(target_path) == False):
        raise FileNotFoundError("The given file path " + str(target_path) + " doesn't exist!")

    if (isdir(destination_path) == False):
        raise ValueError("The given file path " + str(destination_path) + " isn't a directory!")

    if (isdir(target_path) == False):
        raise ValueError("The given file path " + str(target_path) + " isn't a directory!")

    dirpaths = list([target_path])
    filepaths = list([])
    working_directory = wdir()
    
    for root, dirs, files in os.walk(target_path, topdown=True):
        for path in dirs:
            dirpaths.append(str(os.path.join(root, path)))
        for path in files:
            filepaths.append(str(os.path.join(root, path)))

    dirpaths.sort(); filepaths.sort()
    
    chdir(destination_path)

    for dirpath in dirpaths:
        mkdir(dirpath)

    chdir(working_directory)

    for filepath in filepaths:
        copyfile(destination_path, filepath)

# Copy single/multiple files to destination directory.
# The destination_path and *filepaths arguments must be strings.
def copyfile(destination_path, *filepaths):
    working_directory = wdir()

    if (type(destination_path) is not str):
        raise TypeError("The given file path " + str(destination_path) + " isn't a string!")

    if (path_exists(destination_path) == False):
        raise FileNotFoundError("The given file path " + str(destination_path) + " doesn't exist!")

    if (isdir(destination_path) == False):
        raise ValueError("The given file path " + str(destination_path) + " isn't a directory!")

    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (path_exists(filepath) == False):
            raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

        if (isfile(filepath) == False):
            raise ValueError("The given file path " + str(filepath) + " isn't a file!")
  
        data = read(filepath)
        mkfile("b", os.path.join(destination_path, filepath))
        write(os.path.join(destination_path, filepath), data)

# Read the binary from a file and return.
# The filepath argument must be a string.
def read(filepath):
    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isfile(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a file!")

    with open(filepath, "rb") as new_file:
        data = new_file.read(); new_file.close()
        return data

# Write bytes object to a file.
# The filepath argument must be a string.
# The data argument must be a bytes object.
def write(filepath, data):
    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isfile(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a file!")

    if (type(data) is not bytes):
        raise TypeError("The data argument isn't of type bytes!")

    with open(filepath, "wb") as new_file:
        new_file.write(data); new_file.close()
       

# Read and print lines in single/multiple text/binary based files.
# The mode argument must be either strings: "t" (text) or "b" (binary).
# The *filepaths arguments must be strings. 
def rpfile(mode, *filepaths):
    for filepath in filepaths:         
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (type(mode) is not str):
            raise TypeError("The given mode " + str(mode) + " isn't a string!")

        if (path_exists(filepath) == False):
            raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

        if (isfile(filepath) == False):
            raise ValueError("The given file path " + str(filepath) + " isn't a file!")

        if (mode != "t" and mode != "b"):
            raise ValueError("The mode is of the incorrect value! It must either strings: 't' (text) or 'b' (binary)")

        print(" ------------------------------", "\n", "Reading/Printing:", filepath, "\n", "------------------------------")

        with open(filepath, "r" + mode) as new_file:
            for line in new_file:
                print(line, end="")
            
            new_file.close()
            print()

# Delete single/multiple files.
# The *filepaths arguments must be strings.
def rmfile(*filepaths):
    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (path_exists(filepath) == False):
            raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

        if (isfile(filepath) == False):
            raise ValueError("The given file path " + str(filepath) + " isn't a file!")

        os.remove(filepath)

# Make single/multiple text/binary based files.
# The mode argument must be either strings: "t" (text) or "b" (binary).
# The *filepaths arguments must be strings.
def mkfile(mode, *filepaths):
    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")
            
        if (type(mode) is not str):
            raise TypeError("The given mode " + str(mode) + " isn't a string!")

        if (path_exists(filepath) == True):
            raise FileExistsError("The given file path " + str(filepath) + " already exists!")

        if (mode != "t" and mode != "b"):
            raise ValueError("The mode is of the incorrect value! It must either strings: 't' (text) or 'b' (binary)")
        
        new_file = open(filepath, "x" + mode); new_file.close()

# Delete single/multiple empty directories.
# The *filepaths arguments must be strings.
def rmdir(*filepaths):
    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (path_exists(filepath) == False):
            raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

        if (isdir(filepath) == False):
            raise ValueError("The given file path " + str(filepath) + " isn't a directory!")
            
        os.rmdir(filepath)

# Delete single directory along with it's subdirectories and files.
# Use this with caution, as you could delete your entire file system if you're not careful. 
def rmall(dirpath):
    if (type(dirpath) is not str):
        raise TypeError("The given file path " + str(dirpath) + " isn't a string!")

    if (path_exists(dirpath) == False):
        raise FileNotFoundError("The given file path " + str(dirpath) + " doesn't exist!")

    if (isdir(dirpath) == False):
        raise ValueError("The given file path " + str(dirpath) + " isn't a directory!")

    dirpaths = list([])
    filepaths = list([])
    
    for root, dirs, files in os.walk(dirpath, topdown=True):
        for path in dirs:
            dirpaths.append(str(os.path.join(root, path)))
        for path in files:
            filepaths.append(str(os.path.join(root, path)))

    dirpaths.reverse(); filepaths.sort()

    for file in filepaths:
        rmfile(file)

    try:
        for file in dirpaths:
            rmdir(file)
    finally:
        rmdir(dirpath)
 
# Make single/multiple directories.
# The *filepaths arguments must be strings.
def mkdir(*filepaths):
    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (path_exists(filepath) == True):
            raise FileExistsError("The given file path " + str(filepath) + " already exists!")
            
        os.mkdir(filepath)

# Output to terminal the file/subdirectory names of single/multiple argument filepaths.
# Use no arguments for working directory only.
# The *filepaths arguments must be strings. 
def rpdir(*filepaths):
    if (len(filepaths) == 0):
        directory = list(lsdir())
        print(" -------------------", "\n", "Directory - Working", "\n", "-------------------")

        for file in directory:
            print("-", file)        
    else: 
        for filepath in filepaths:
            if (type(filepath) is not str):
                raise TypeError("The given file path " + str(filepath) + " isn't a string!")

            if (path_exists(filepath) == False):
                raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

            if (isdir(filepath) == False):
                raise ValueError("The given file path " + str(filepath) + " isn't a directory!")
                
            directory = list(lsdir(filepath))
            print(" -------------------", "\n", "Directory - ", filepath.strip("../"), "\n", "-------------------")

            for file in directory:
                print("-", file)

# Return a list with file/subdirectory names of the single/multiple argument filepaths.
# If there are no arguments used in *filepaths, a list of the contents within the working directory is returned.
# If there is only one argument used in *filepaths, a list of the contents of only that directory is returned. 
# Using the method with two or more arguments in *filepaths will return a list of lists with each list containing the file/subdirectory names of that filepath argument.
def lsdir(*filepaths):
    if (len(filepaths) == 0):
        return list(os.listdir())
    elif (len(filepaths) == 1):
        if (type(filepaths[0]) is not str):
            raise TypeError("The given file path " + str(filepaths[0]) + " isn't a string!")

        if (path_exists(filepaths[0]) == False):
            raise FileNotFoundError("The given file path " + str(filepaths[0]) + " doesn't exist!")

        if (isdir(filepaths[0]) == False):
            raise ValueError("The given file path " + str(filepaths[0]) + " isn't a directory!")

        return list(os.listdir(filepaths[0]))
    else:
        finalList = list()

        for filepath in filepaths:
            fileList = list()

            if (type(filepath) is not str):
                raise TypeError("The given file path " + str(filepath) + " isn't a string!")

            if (path_exists(filepath) == False):
                raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

            if (isdir(filepath) == False):
                raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

            for file in os.listdir(filepath):
                fileList.append(file)

            finalList.append(fileList)

        return finalList

# Change current working directory.
# The filepath argument must be a string. 
def chdir(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False and filepath != ".."):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isdir(filepath) == False and filepath != ".."):
        raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

    os.chdir(filepath)

# Return string of the path of the current working directory. 
def wdir():
    return str(os.getcwd())

# Print working directory to terminal. 
def pwdir():
    print("-", wdir())

# Return a list from lines in single/multiple text based files.
# The *filepaths arguments must be strings.
def mklist(*filepaths):
    finalList = list([])
    tempList = list([])
        
    for filepath in filepaths:
        if (type(filepath) is not str):
            raise TypeError("The given file path " + str(filepath) + " isn't a string!")

        if (path_exists(filepath) == False):
            raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

        if (isfile(filepath) == False):
            raise ValueError("The given file path " + str(filepath) + " isn't a file!")    

        with open(filepath, "rt") as new_file:
            tempList = new_file.readlines()
            new_file.close()

            for element in tempList:
                element = element.strip("\n")
                finalList.append(element)
               
    return finalList
 
# Write singular strings or lists of strings in sequence to lines in a text based file.
# The filepath argument must be a string. 
# The lines in the file are overwritten by the lines argument values.       
def writelines(filepath, *lines):
    j = int(0)
    i = int(0)

    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isfile(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a file!")

    with open(filepath, "wt") as new_file:
        for line in lines:      
            if (type(line) is not str and type(line) is not list):
                raise ValueError("The line argument must be a string or a list of strings!")

            if (type(line) is str):
                if (line == ""):
                    pass
                elif (j == 0 and line != ""):
                    new_file.write(line)
                else:
                    new_file.write("\n")
                    new_file.write(line)
                        
                j += 1
                    
            if (type(line) is list):
                if (len(line) == 0):
                    pass
                elif (j == 0 and len(line) != 0):

                    while (i < len(line) - 1):
                        if (type(line[i]) is not str):
                            raise TypeError("A value in the list argument isn't a string!")

                        new_file.write(line[i])
                        new_file.write("\n")
                        i += 1

                    new_file.write(line[i])
                    i = int(0)
                else:
                    new_file.write("\n")

                    while (i < len(line) - 1):
                        if (type(line[i]) is not str):
                            raise TypeError("A value in the list argument isn't a string!")

                        new_file.write(line[i])
                        new_file.write("\n")
                        i += 1

                    new_file.write(line[i])
                    i = int(0)

                j += 1
                
        new_file.close()

# Append singular strings or lists of strings in sequence to lines at the end of a text based file. 
# The filepath argument must be a string.
def appendlines(filepath, *lines):
    i = int(0)

    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")
     
    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isfile(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a file!")

    with open(filepath, "at") as new_file:
        for line in lines:      
            if (type(line) is not str and type(line) is not list):
                raise ValueError("The line argument must be a string or a list of strings!")

            if (type(line) is str):        
                new_file.write("\n")
                new_file.write(line)
                        
            if (type(line) is list):
                if (len(line) == 0):
                    pass
                else: 
                    new_file.write("\n")
                    while (i < len(line) - 1):
                        new_file.write (str(line[i]))
                        new_file.write("\n")
                        i += 1

                    new_file.write(line[i])
                    i = int(0)             
                                         
        new_file.close()

# Return boolean value (True or False) to check if a single file path exists.
# The filepath argument must be a string.
def path_exists(filepath):  
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (os.path.exists(filepath) == True):
        return True
    elif (os.path.exists(filepath) == False):
        return False

# Return boolean value (True or False) to check if filepath argument is a file.
# The filepath argument must be a string.
def isfile(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (os.path.isfile(filepath) == True):
        return True
    elif (os.path.isfile(filepath) == False):
        return False

# Return boolean value (True or False) to check if filepath argument is a directory.
# The filepath argument must be a string.
def isdir(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (os.path.isdir(filepath) == True):
        return True
    elif (os.path.isdir(filepath) == False):
        return False

# Return boolean value (True or False) to check if filepath argument is a tar archive.
# The filepath argument must be a string.
def istar(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (tarfile.is_tarfile(filepath) == True):
        return True
    elif (tarfile.is_tarfile(filepath) == False):
        return False

# Return boolean value (True or False) to check if filepath argument is a zip archive.
# The filepath argument must be a string.
def iszip(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (zipfile.is_zipfile(filepath) == True):
        return True
    elif (zipfile.is_zipfile(filepath) == False):
        return False

'''
    Rename single/multiple files or directories.
    
    current_filepath represents the file path's name being changed.
    desired_filepath represents the file path's new intended name. 
    
    current_filepath and desired_filepath can either be:
        Two strings
        Two lists of equal length consisting of strings
'''
def ren(current_filepath, desired_filepath):
    i = int(0)

    if (type(current_filepath) is str and type(desired_filepath) is not str):
        raise TypeError("The current_filepath argument is a string, the desired_filepath argument must be a string!")

    if (type(current_filepath) is not str and type(desired_filepath) is str):
        raise TypeError("The desired_filepath argument is a string, the current_filepath argument must be a string!")

    if (type(current_filepath) is list and type(desired_filepath) is not list):
        raise TypeError("The current_filepath argument is a list, the desired_filepath argument must be a list!")

    if (type(current_filepath) is not list and type(desired_filepath) is list):
        raise TypeError("The desired_filepath argument is a list, the current_filepath argument must be a list!")
        
    if (type(current_filepath) is str and type(desired_filepath) is str):
        if (path_exists(current_filepath) == False):
            raise FileNotFoundError("The given file path " + str(current_filepath) + " doesn't exist!")

        if (path_exists(desired_filepath) == True):
            raise FileExistsError("The given file path " + str(desired_filepath) + " already exists!")

        os.rename(current_filepath, desired_filepath)

    if (type(current_filepath) is list and type(desired_filepath) is list):
        if (len(current_filepath) != len(desired_filepath)):
            raise ValueError("The length of the list current_filepath is not equal to the length of the list desired_filepath!")

        while(i < len(current_filepath) and i < len(desired_filepath)):
            if (type(current_filepath[i]) is not str):
                raise TypeError("The given file path " + str(current_filepath[i]) + " isn't a string!")

            if (type(desired_filepath[i]) is not str):
                raise TypeError("The given file path " + str(desired_filepath[i]) + " isn't a string!")

            if (path_exists(current_filepath[i]) == False):
                raise FileNotFoundError("The given file path " + str(current_filepath[i]) + " doesn't exist!")

            if (path_exists(desired_filepath[i]) == True):
                raise FileExistsError("The given file path " + str(desired_filepath[i]) + " already exists!")

            os.rename(current_filepath[i], desired_filepath[i])
            i += 1

# Create a tar archive with gzip compression & .gz extension.
def tar_wrap(filepath):
    arcpath = str(filepath + ".tar.gz")
    dirpaths = list([])
    filepaths = list([])

    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isdir(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

    if (path_exists(arcpath) == True):
        raise FileExistsError("The archive" + arcpath + "already exists!")

    for root, dirs, files in os.walk(filepath, topdown=True):
        for path in dirs:
            dirpaths.append(str(os.path.join(root, path)))
        for path in files:
            filepaths.append(str(os.path.join(root, path)))

    dirpaths.sort(); filepaths.sort()

    with tarfile.open(name=arcpath, mode='w:gz') as archive:
        for file in dirpaths:
            archive.add(file, recursive=False)

        for file in filepaths:
            archive.add(file, recursive=False)

        archive.close()

# Extract a tar gzip archive contents to working directory.
def tar_extract(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (istar(filepath) == False or str(filepath).endswith(".tar.gz") == False):
        raise FileNotFoundError("The given file path isn't a .tar.gz archive!")

    if (path_exists(str(filepath).strip(".tar.gz")) == True):
        raise FileExistsError("The extraction directory already exists!")

    with tarfile.open(name=filepath, mode='r:gz') as archive:
        archive.extractall()
        archive.close()

# Create a zip archive with DEFLATE compression.
def zip_wrap(filepath):
    arcpath = str(filepath + ".zip")
    dirpaths = list([])
    filepaths = list([])

    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isdir(filepath) == False):
        raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

    if (path_exists(arcpath) == True):
        raise FileExistsError("The archive" + arcpath + "already exists!")

    for root, dirs, files in os.walk(filepath, topdown=True):
        for path in dirs:
            dirpaths.append(str(os.path.join(root, path)))
        for path in files:
            filepaths.append(str(os.path.join(root, path)))

    dirpaths.sort(); filepaths.sort()

    with zipfile.ZipFile(arcpath, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        for file in dirpaths:
            archive.write(file)

        for file in filepaths:
            archive.write(file)

        archive.close()

# Extract a zip archive contents to working directory.
def zip_extract(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (iszip(filepath) == False or str(filepath).endswith(".zip") == False):
        raise FileNotFoundError("The given file path isn't a .zip archive!")

    if (path_exists(str(filepath).strip(".tar.gz")) == True):
        raise FileExistsError("The extraction directory already exists!")

    with zipfile.ZipFile(filepath, mode="r", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.extractall()
        archive.close()

# Count and return the number of files within a directory.
def filecount(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isdir(filepath) == False and filepath != ".."):
        raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

    working_directory = wdir()
    count = int(0)
    contents = list(lsdir(filepath))
    chdir(filepath)

    for item in contents:
        if (isfile(item) == True):
            count += 1
    
    chdir(working_directory)

    return count

# Count and return the number of directories within a directory.
def dircount(filepath):
    if (type(filepath) is not str):
        raise TypeError("The given file path " + str(filepath) + " isn't a string!")

    if (path_exists(filepath) == False):
        raise FileNotFoundError("The given file path " + str(filepath) + " doesn't exist!")

    if (isdir(filepath) == False and filepath != ".."):
        raise ValueError("The given file path " + str(filepath) + " isn't a directory!")

    working_directory = wdir()
    count = int(0)
    contents = list(lsdir(filepath))
    chdir(filepath)

    for item in contents:
        if (isdir(item) == True):
            count += 1
    
    chdir(working_directory)

    return count