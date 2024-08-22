# Depends: stringr package

# Read a single OceanView ASCII file
read.OceanView = function(file, skip=14){

  # Get the irradiance data
  measures = read.table(file, skip=skip)
  
  # Parse header metadata
  header = read.table(file, sep='\n', nrow=skip)

  ## Metadata extract variables
  folder_depth = 1 + stringr::str_count(file, pattern='/')
  fn = stringr::str_split_i(file, pattern = '/', i=folder_depth) #This allows us to isolate the filename for the filename column in the dataframe, without the full filepath.
  
  time = stringr::str_extract(header[2,], pattern='\\d{2}\\:\\d{2}\\:\\d{2}')
  integration = stringr::str_extract(header[6,], pattern='(\\d+\\.\\d+E-\\d*)|(\\d+\\.\\d+)|(\\d+)') #Integration time setting when running the spectrometer
  scans = stringr::str_extract(header[7,], pattern='\\d+') #How many scans to average
  
  rm(folder_depth)

  # Format metadata - use nrows of measures to make big vectors of the metadata
  filename = rep(fn, nrow(measures))
  time_vec = rep(time, nrow(measures))
  int_vec = rep(integration, nrow(measures))
  scans_vec = rep(scans, nrow(measures))
  
  metadata = cbind(filename, time_vec, int_vec, scans_vec)

  # Make dataframe
  combined = data.frame(cbind(metadata, measures))
  colnames(combined) = c('filename', 'time', 'integration_time', 'scans', 'wavelength', 'irradiance')
  return(combined)
}

#---
# Read a folder of OceanView ASCII files

read_many.OceanView = function(filepath){
  filenames = list.files(filepath, pattern='.txt') #Make a list of filenames

  # Go through each file and extract the data into a list
  measurements = lapply(filenames, function(f){
    print(f) # Reassures user when importing lots of files
    full_path = paste(filepath, f, sep='/') #Reconstruct the full path to pass into read.OceanView
    data = read.OceanView(full_path)
    return(data)
  })

  # Formatting
  measurements = as.data.frame(do.call(rbind, measurements)) #Converts from list to dataframe
  print(summary(measurements)) #Gives user info about their data
  return(measurements)
}

#---
# Match the timestamp with the event number from the light schedule

event_nos_timestamp = function(){ 
  # Wil write at a later stage. Adapt from HelioProgramming test 20
}
