# Extrapolating data to Indian coasts: Parallel processing
library(ncdf4)

is_in_polygon = function(ch_lat, ch_lon)
{
  library(sf)
  map = read_sf("Data/2025/Shp files/SA_countries.shp")
  
  pnts <- data.frame(
    "x" = ch_lon,
    "y" = ch_lat)
  
  pnts_sf <- sf::st_as_sf(pnts, coords = c("x", "y"), crs = st_crs(map))
  xx=st_intersects(pnts_sf, map)
  yy=ifelse(unlist(lapply(xx, length))==0,TRUE,FALSE)
  return(yy)
}

#Function to fill empty cells
fill_cell=function(vals, empty_cells, lat, lon, time, lat_lb, lat_ub, lon_lb, lon_ub)#water_v
{
  for(la in 1:length(lat_lb:lat_ub))
  {
    if(length(empty_cells[[la]])!=0)  
    {
      for(c in 1:length(empty_cells[[la]]))
      {
        cell_lon=which(lon==empty_cells[[la]][c])
        cell_lat=la+lat_lb-1
        
        # if(cell_lat==lat[1])
        #   cu4=cu5=cu6=NA
        # if(cell_lat==lat[length(lat)])
        #   cu2=cu1=cu8=NA
        # if(cell_lon==lon[1])
        #   cu2=cu3=cu4=NA
        # if(cell_lon==lon[length(lon)])
        #   cu6=cu7=cu8=NA
        
        cu1=vals[cell_lon, cell_lat+1]
        cu2=vals[cell_lon-1, cell_lat+1]
        cu3=vals[cell_lon-1, cell_lat]
        cu4=vals[cell_lon-1, cell_lat-1]
        cu5=vals[cell_lon, cell_lat-1]
        cu6=vals[cell_lon+1, cell_lat-1]
        cu7=vals[cell_lon+1, cell_lat]
        cu8=vals[cell_lon+1, cell_lat+1]
        cell_val_u=mean(c(cu1,cu2,cu3,cu4,cu5,cu6,cu7,cu8), na.rm=T)
        vals[cell_lon, cell_lat]=round(cell_val_u,3)
      }
    }
  }
  return(vals)
}

wd="D:/Quantifying cetacean bycatch/HYCOM/NC_files/"
files=list.files(wd, ".nc")
index = index = grep("20240901T0000Z", files)
files=files[1:(index-1)]
#files=files[index:length(files)]

#Finding out empty cells
filename=files[1]
ncdata=nc_open(paste0(wd, filename), write=F)
# print(ncdata)
# nc_close(ncdata)
lon=round(ncvar_get(ncdata,"lon"),4)
#lon=lon[-c(1,length(lon))]
lat=round(ncvar_get(ncdata,"lat"),4)
#lat=lat[-c(1,length(lat))]
water_u=ncvar_get(ncdata,"water_u")
time=ncvar_get(ncdata,"time")
nc_close(ncdata)

#extent where extrapolation is needed
lat_lb=which.min(abs(lat-5))
lat_ub=which.min(abs(lat-25))
lon_lb=which.min(abs(lon-68))
lon_ub=which.min(abs(lon-98))

#Finding out empty cells; Skip edge cells as they require a modified algorithm for filling

#wateru=water_u[,,1]#2:length(lon),2:length(lat) #for multidimension (time) data
wateru=water_u#2:length(lon),2:length(lat)
empty_cells=list()
stm = Sys.time()
for(la in 1:length(lat_lb:lat_ub))
{
  na_lon=lon[which(is.na(wateru[,la+lat_lb-1]))] # lon, lat, time la+lat_lb-1
  na_lon=na_lon[na_lon>=lon[lon_lb] & na_lon<=lon[lon_ub]]
  if(length(na_lon)!=0)
  {
    na_lon=na_lon[sapply(lat[la+lat_lb-1], na_lon, FUN=is_in_polygon)] # get empty cells that are in water
    if(lon[1] %in% na_lon == T)# remove edge longitude cells
      na_lon=na_lon[-which(na_lon==lon[1])]
    if(lon[length(lon)] %in% na_lon == T)# remove edge longitude cells
      na_lon=na_lon[-which(na_lon==lon[length(lon)])]
  }
  empty_cells[[la]]=na_lon
}
Sys.time()

#Discarding empty cells in first and last latitude rows for easy computation
empty_cells[[1]]=numeric(0)
empty_cells[[length(lat)]]=numeric(0)

for(i in 1:length(empty_cells))
{
  if(lon[1] %in% empty_cells[[i]] == T)# remove edge longitude cells
    empty_cells[[i]]=empty_cells[[i]][-which(empty_cells[[i]]==lon[1])]
  if(lon[length(lon)] %in% na_lon == T)# remove edge longitude cells
    empty_cells[[i]]=empty_cells[[i]][-which(empty_cells[[i]]==lon[length(lon)])]
}
Sys.time() - stm

#Fill empty cells with mean values
stm = Sys.time()
for(h in 1:length(files))
{
  #Finding out empty cells
  filename=files[h]
  ncdata=nc_open(paste0(wd, filename), write=F)
  # print(ncdata)
  # nc_close(ncdata)
  lon=round(ncvar_get(ncdata,"lon"),4)
  #lon=lon[-c(1,length(lon))]
  lat=round(ncvar_get(ncdata,"lat"),4)
  #lat=lat[-c(1,length(lat))]
  water_u=ncvar_get(ncdata,"water_u")
  time=ncvar_get(ncdata,"time")
  nc_close(ncdata)
 
  #Formatting data for entry
  # For multidimension data
  # vals_u=list()
  # for(i in 1:dim(water_u)[3])
  #   vals_u[[i]]=(as.matrix(water_u[,,i]))
  # 
  # vals_v=list()
  # for(i in 1:dim(water_v)[3])
  #   vals_v[[i]]=(as.matrix(water_v[,,i]))
  
  #For single dimension data
  vals_u=list(water_u)
  vals_v=list(water_v)
  
  # library(parallel)
  # cl <- makeCluster(12)
  # #clusterExport(cl, "is_in_polygon")
  # u=parLapply(cl, vals_u, fill_cell, empty_cells, lat, lon, time, lat_lb, lat_ub, lon_lb, lon_ub)
  # v=parLapply(cl, vals_v, fill_cell, empty_cells, lat, lon, time, lat_lb, lat_ub, lon_lb, lon_ub)
  # stopCluster(cl)
  
  u=lapply(vals_u, FUN=fill_cell, empty_cells, lat, lon, time, lat_lb, lat_ub, lon_lb, lon_ub)
  v=lapply(vals_v, FUN=fill_cell, empty_cells, lat, lon, time, lat_lb, lat_ub, lon_lb, lon_ub)
  
  #Formatting data for re-entry
  # final_vals_u=array(numeric(),c(nrow(u[[1]]),(ncol(u[[1]])),length(u)))
  # for(i in 1:length(u))
  #   final_vals_u[,,i]=u[[i]]
  # final_vals_v=array(numeric(),c(nrow(v[[1]]),(ncol(v[[1]])),length(v)))
  # for(i in 1:length(v))
  #   final_vals_v[,,i]=v[[i]]
  
  final_vals_u=u[[1]]
  final_vals_v=v[[1]]
  
  # write data to a new file because the original can't be appended for some reason
  
  # Define the dimensions and variables
  dimX = ncdim_def("lon", "degrees", lon)
  dimY = ncdim_def("lat", "degrees", lat)
  dimT = ncdim_def("time", "hours since 2000-01-01", time)
  
  
  mv=-30000 # missing value
  water_u_val=ncvar_def('water_u', 'm/s', list(dimX,dimY,dimT), mv, prec="float")
  water_v_val=ncvar_def('water_v', 'm/s', list(dimX,dimY,dimT), mv, prec="float")
  
  # Create the NetCDF file
  my_new_nc = nc_create(paste0(wd,"updated/udt_",filename),  list(water_u_val, water_v_val))
  
  #print(my_new_nc)
  # Write data to the NetCDF file
  ncvar_put(my_new_nc, "water_u", final_vals_u)
  ncvar_put(my_new_nc, "water_v", final_vals_v)
  
  # Close the new file to finish writing
  nc_close(my_new_nc)
}
Sys.time() - stm