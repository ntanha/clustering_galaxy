# routines for reading headers and data blocks from Gadget snapshot files
# usage e.g.:
#
# import readsnap as rs
# header = rs.snapshot_header("snap_063.0") # reads snapshot header
# print header.massarr
# mass = rs.read_block("snap_063","MASS",parttype=5) # reads mass for particles of type 5, using block names should work for both format 1 and 2 snapshots
# print "mass for", mass.size, "particles read"
# print mass[0:10]
#
# before using read_block, make sure that the description (and order if using format 1 snapshot files) of the data blocks
# is correct for your configuration of Gadget 
#
# for mutliple file snapshots give e.g. the filename "snap_063" rather than "snap_063.0" to read_block
# for snapshot_header the file number should be included, e.g."snap_063.0", as the headers of the files differ
#
# the returned data block is ordered by particle species even when read from a multiple file snapshot

import numpy as np
import os
import sys
import math
  
# ----- class for snapshot header ----- 

class snapshot:
  def __init__(self, filename=None, **kwargs):
    if filename != None:
      self.read(filename, **kwargs)
    
  def read(self, filename, export=False):
    self.head = snapshot_header(filename)
    self.ntot = sum(self.head.npart)
    self.pos  = read_block(filename, "POS ")  
    self.vel  = read_block(filename, "VEL ")  
    self.ids  = read_block(filename, "ID  ")  
    self.mass = read_block(filename, "MASS")  
    self.u    = read_block(filename, "U   ")
    self.rho  = read_block(filename, "RHO ")
    if self.head.cooling:
      self.ne   = read_block(filename, "NE  ")
      self.nh   = read_block(filename, "NH  ")
    self.hsml = read_block(filename, "HSML")
    if self.head.sfr:
      self.sfr  = read_block(filename, "SFR ")
      self.m_imf= read_block(filename, "MIMF")
    if self.head.age:
      self.age  = read_block(filename, "AGE ")
    if self.head.metals:
      self.let = read_block(filename, "LET ")
      self.imass = read_block(filename, "INIM")
      self.metgas = read_block(filename, "Z   ", parttype=0, csformat = 1)
      self.metstars = read_block(filename, "Z   ", parttype=4, csformat = 1)
      self.temp = read_block(filename, "CSTE", csformat = 1)
      self.pot = read_block(filename, "POT ", csformat = 1)
    

    parttypes = ['gas', 'halo', 'disk', 'bulge', 'stars', 'bh']
    startind = 0
    for i in xrange(len(parttypes)):
      endind = startind + self.head.npart[i]
      varname = ''.join(["n", parttypes[i]])
      self.__dict__[varname] = self.head.npart[i]
      varname = ''.join(["x", parttypes[i]])
      self.__dict__[varname] = self.pos[startind:endind, 0]
      varname = ''.join(["y", parttypes[i]])
      self.__dict__[varname] = self.pos[startind:endind, 1]
      varname = ''.join(["z", parttypes[i]])
      self.__dict__[varname] = self.pos[startind:endind, 2]
      varname = ''.join(["vx", parttypes[i]])
      self.__dict__[varname] = self.vel[startind:endind, 0]
      varname = ''.join(["vy", parttypes[i]])
      self.__dict__[varname] = self.vel[startind:endind, 1]
      varname = ''.join(["vz", parttypes[i]])
      self.__dict__[varname] = self.vel[startind:endind, 2]
      varname = ''.join(["m", parttypes[i]])
      self.__dict__[varname] = self.mass[startind:endind]
      varname = ''.join(["i", parttypes[i]])
      self.__dict__[varname] = self.ids[startind:endind]
      startind = endind

    if export:
      self.export()

  def calcTemp(self):
      gamma_minus_1 = 5.0/3 -1
      boltzmann = 1.3806e-16
      protonm= 1.6726e-24
      hydr_frac= 0.76
      yhelium = ( 1 - hydr_frac ) / ( 4 * hydr_frac )
      mass_in_g= 1.989e43
      length_cm=3.085678e21
      vel_in_cm_per_s=1e5
      time_in_s= length_cm / vel_in_cm_per_s
      energy_cgs = mass_in_g * (length_cm**2) / (time_in_s**2);
      mu = (1 + 4 * yhelium) / (1 + yhelium + self.ne);          
      temp = gamma_minus_1 / boltzmann * self.u * protonm * mu;
      temp *= energy_cgs / mass_in_g;
      self.temp = temp

  def export(self):
    for k in self.__dict__:
      sys.modules['__builtin__'].__dict__[k] = self.__dict__[k]
    

class snapshot_header:
  def __init__(self, filename):
    if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()
      
    self.filename = filename  
    f = open(filename,'rb')    
    blocksize = np.fromfile(f,dtype=np.int32,count=1)
    if blocksize[0] == 8:
      swap = 0
      format = 2
    elif blocksize[0] == 256:
      swap = 0
      format = 1  
    else:
      blocksize.byteswap(True)
      if blocksize[0] == 8:
        swap = 1
        format = 2
      elif blocksize[0] == 256:
        swap = 1
        format = 1
      else:
        print "incorrect file format encountered when reading header of", filename
        sys.exit()
    
    self.format = format
    self.swap = swap
    
    if format==2:
      f.seek(16, os.SEEK_CUR)
    
    self.npart = np.fromfile(f,dtype=np.int32,count=6)
    self.massarr = np.fromfile(f,dtype=np.float64,count=6)
    self.time = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.redshift = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.sfr = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.feedback = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.nall = np.fromfile(f,dtype=np.int32,count=6)
    self.cooling = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.filenum = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.boxsize = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_m = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_l = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.hubble = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.age =  (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.metals =  (np.fromfile(f,dtype=np.int32,count=1))[0]
    
    if swap:
      self.npart.byteswap(True)
      self.massarr.byteswap(True)
      self.time = self.time.byteswap()
      self.redshift = self.redshift.byteswap()
      self.sfr = self.sfr.byteswap()
      self.feedback = self.feedback.byteswap()
      self.nall.byteswap(True)
      self.cooling = self.cooling.byteswap()
      self.filenum = self.filenum.byteswap()
      self.boxsize = self.boxsize.byteswap()
      self.omega_m = self.omega_m.byteswap()
      self.omega_l = self.omega_l.byteswap()
      self.hubble = self.hubble.byteswap()
      self.age =  self.age.byteswap()
      self.metals =  self.metals.byteswap()
     
    f.close()
 
# ----- find offset and size of data block ----- 

def find_block(filename, format, swap, block, block_num, only_list_blocks=False):
  if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()
            
  f = open(filename,'rb')
  f.seek(0, os.SEEK_END)
  filesize = f.tell()
  f.seek(0, os.SEEK_SET)
  
  found = False
  curblock_num = 1
  while ((not found) and (f.tell()<filesize)):
    if format==2:
      f.seek(4, os.SEEK_CUR)
      curblock = f.read(4)
      if (block == curblock):
        found = True
      f.seek(8, os.SEEK_CUR)  
    else:
      if curblock_num==block_num:
        found = True
        
    curblocksize = (np.fromfile(f,dtype=np.int32,count=1))[0]
    if swap:
      curblocksize = curblocksize.byteswap()
    
    # - print some debug info about found data blocks -
    #if format==2:
    #  print curblock, curblock_num, curblocksize
    #else:
    #  print curblock_num, curblocksize
    
    if only_list_blocks:
      print curblock_num,curblock,f.tell(),curblocksize
      found = False
    
    if found:
      blocksize = curblocksize
      offset = f.tell()
    else:
      f.seek(curblocksize, os.SEEK_CUR)
      blocksize_check = (np.fromfile(f,dtype=np.int32,count=1))[0]
      if swap: blocksize_check = blocksize_check.byteswap()
      if (curblocksize != blocksize_check):
        print "something wrong"
        sys.exit()
      curblock_num += 1
      
  f.close()
      
  if ((not found) and (not only_list_blocks)):
    print "Error: block not found"
    sys.exit()
    
  if (not only_list_blocks):
    return offset,blocksize
 
# ----- read data block -----
 
def read_block(filename, block, parttype=-1, physical_velocities=True, arepo=0, no_masses=False, verbose=False, csformat=0):
  if (verbose):
	  print "reading block", block
  
  blockadd=0
  blocksub=0
  
  if arepo==0:
    if (verbose):	
	    print "Gadget format"
    blockadd=0
  if arepo==1:
    if (verbose):	
	    print "Arepo format"
    blockadd=1	
  if arepo==2:
    if (verbose):
	   print "Arepo extended format"
    blockadd=4	
  if no_masses==True:
    if (verbose):	
	    print "No mass block present"    
    blocksub=1
		 
  if parttype not in [-1,0,1,2,3,4,5]:
    print "wrong parttype given"
    sys.exit()
  
  if os.path.exists(filename):
    curfilename = filename
  elif os.path.exists(filename+".0"):
    curfilename = filename+".0"
  else:
    print "file not found:", filename
    print "and:", curfilename
    sys.exit()
  
  head = snapshot_header(curfilename)
  format = head.format
  swap = head.swap
  npart = head.npart
  massarr = head.massarr
  nall = head.nall
  filenum = head.filenum
  redshift = head.redshift
  time = head.time
  del head
  
  # - description of data blocks -
  # add or change blocks as needed for your Gadget version
  data_for_type = np.zeros(6,bool) # should be set to "True" below for the species for which data is stored in the data block
  dt = np.float32 # data type of the data in the block
  if block=="POS ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 2
  elif block=="VEL ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 3
  elif block=="ID  ":
    data_for_type[:] = True
    dt = np.uint32
    block_num = 4
  elif block=="MASS":
    data_for_type[np.where(massarr==0)] = True
    block_num = 5
    if parttype>=0 and massarr[parttype]>0:   
      if (verbose):	
	      print "filling masses according to massarr"   
      return np.ones(nall[parttype],dtype=dt)*massarr[parttype]
  elif block=="U   ":
    data_for_type[0] = True
    block_num = 6-blocksub
  elif block=="RHO ":
    data_for_type[0] = True
    block_num = 7-blocksub
  elif block=="VOL ":
    data_for_type[0] = True
    block_num = 8-blocksub 
  elif block=="CMCE":
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 9-blocksub 
  elif block=="AREA":
    data_for_type[0] = True
    block_num = 10-blocksub
  elif block=="NFAC":
    data_for_type[0] = True
    dt = np.dtype(np.int32)	
    block_num = 11-blocksub
  elif block=="NE  ":
    data_for_type[0] = True
    block_num = 8+blockadd-blocksub
  elif block=="NH  ":
    data_for_type[0] = True
    block_num = 9+blockadd-blocksub
  elif block=="HSML":
    data_for_type[0] = True
    block_num = 10+blockadd-blocksub
  elif block=="SFR ":
    data_for_type[0] = True
    block_num = 11+blockadd-blocksub
  elif block=="AGE ":
    data_for_type[4] = True
    block_num = 12+blockadd-blocksub
  elif block=="LET ":
    data_for_type[4] = True
    block_num = 13+blockadd-blocksub
  elif block=="INIM":
    data_for_type[4] = True
    block_num = 14+blockadd-blocksub
  elif block=="Z   ":
    data_for_type[0] = True
    data_for_type[4] = True
    block_num = 13+blockadd-blocksub
    if csformat:
      dt = np.dtype((np.float32,12))
      block_num = 15+blockadd-blocksub
  elif block=="POT ":
    data_for_type[:] = True
    block_num = 13+blockadd-blocksub
    if csformat:
      block_num = 16+blockadd-blocksub
  elif block=="CSTE":
    data_for_type[0] = True
    block_num = 17+blockadd-blocksub
#
  elif block=="CSHS":
    data_for_type[4] = True
    block_num = 18+blockadd-blocksub    
  elif block=="COLN":
    dt = np.dtype((np.float32,12))
    data_for_type[0] = True
    block_num = 19+blockadd-blocksub  
  elif block=="CHEM":
    dt = np.dtype((np.float32,3))
    data_for_type[0] = True
    block_num = 20+blockadd-blocksub  
  elif block=="GAMM":
    data_for_type[0] = True
    block_num = 21+blockadd-blocksub      
  elif block=="CHET":
    data_for_type[0] = True
    block_num = 22+blockadd-blocksub     
  elif block=="DUST":
    data_for_type[0] = True
    block_num = 23+blockadd-blocksub         
  elif block=="COOL":    
    data_for_type[0] = True
    block_num = 24+blockadd-blocksub     
  elif block=="CSSI":
    data_for_type[0] = True
    block_num = 25+blockadd-blocksub
  elif block=="G0  ":
    data_for_type[0] = True
    block_num = 26+blockadd-blocksub
  elif block=="SHH2":
    data_for_type[0] = True
    block_num = 27+blockadd-blocksub
  elif block=="SHDU":
    data_for_type[0] = True
    block_num = 28+blockadd-blocksub
  elif block=="CHC1":
    dt = np.dtype((np.float32,28))
    data_for_type[0] = True
    block_num = 29+blockadd-blocksub
  elif block=="CHC2":
    dt = np.dtype((np.float32,6))
    data_for_type[0] = True
    block_num = 30+blockadd-blocksub
  elif block=="ENDT":
    data_for_type[0] = True
    block_num = 31+blockadd-blocksub

  elif block=="SFFF":
    data_for_type[0] = True
    block_num = 32+blockadd-blocksub
  elif block=="UVIS":
    data_for_type[0] = True
    block_num = 33+blockadd-blocksub
  elif block=="UDIF":
    data_for_type[0] = True
    block_num = 34+blockadd-blocksub

  elif block=="BHMA":
    data_for_type[5] = True
    block_num = 14+blockadd-blocksub
  elif block=="BHMD":
    data_for_type[5] = True
    block_num = 15+blockadd-blocksub
  elif block=="COOR":
    data_for_type[0] = True
    block_num = -1
  elif block=="MIMF":
    dt = np.dtype((np.float32,10))
    data_for_type[4] = True
    block_num = 35+blockadd-blocksub
  else:
    print "Sorry! Block type", block, "not known!"
    sys.exit()
  # - end of block description -

  if (block_num < 0 and format==1):
    print "Sorry! Block number of", block, "not known! Unable to read this block from format 1 file!"
    sys.exit() 
    
  actual_data_for_type = np.copy(data_for_type)  
  if parttype >= 0:
    actual_data_for_type[:] = False
    actual_data_for_type[parttype] = True
    if data_for_type[parttype]==False:
      print "Error: no data for specified particle type", parttype, "in the block", block   
      sys.exit()
  elif block=="MASS":
    actual_data_for_type[:] = True  
    
  allpartnum = np.int64(0)
  species_offset = np.zeros(6,np.int64)
  for j in range(6):
    species_offset[j] = allpartnum
    if actual_data_for_type[j]:
      allpartnum += nall[j]
    
  for i in range(filenum): # main loop over files
    if filenum>1:
      curfilename = filename+"."+str(i)
      
    if i>0:
      head = snapshot_header(curfilename)
      npart = head.npart  
      del head
      
    curpartnum = np.int32(0)
    cur_species_offset = np.zeros(6,np.int64)
    for j in range(6):
      cur_species_offset[j] = curpartnum
      if data_for_type[j]:
        curpartnum += npart[j]
    
    if parttype>=0:
      actual_curpartnum = npart[parttype]      
      add_offset = cur_species_offset[parttype] 
    else:
      actual_curpartnum = curpartnum
      add_offset = np.int32(0)
      
    offset,blocksize = find_block(curfilename,format,swap,block,block_num)
    
    if i==0: # fix data type for ID if long IDs are used
      if block=="ID  ":
        if blocksize == np.dtype(dt).itemsize*curpartnum * 2:
          dt = np.uint64 
        
    if np.dtype(dt).itemsize*curpartnum != blocksize:
      print "something wrong with blocksize! expected =",np.dtype(dt).itemsize*curpartnum,"actual =",blocksize, " in block ", block
      sys.exit()
    
    f = open(curfilename,'rb')
    f.seek(offset + add_offset*np.dtype(dt).itemsize, os.SEEK_CUR)  
    curdat = np.fromfile(f,dtype=dt,count=actual_curpartnum) # read data
    f.close()  
    if swap:
      curdat.byteswap(True)  
      
    if i==0:
      data = np.empty(allpartnum,dt)
    
    for j in range(6):
      if actual_data_for_type[j]:
        if block=="MASS" and massarr[j]>0: # add mass block for particles for which the mass is specified in the snapshot header
          data[species_offset[j]:species_offset[j]+npart[j]] = massarr[j]
        else:
          if parttype>=0:
            data[species_offset[j]:species_offset[j]+npart[j]] = curdat
          else:
            data[species_offset[j]:species_offset[j]+npart[j]] = curdat[cur_species_offset[j]:cur_species_offset[j]+npart[j]]
        species_offset[j] += npart[j]

    del curdat

  if physical_velocities and block=="VEL " and redshift!=0:
    data *= math.sqrt(time)

  return data
  
# ----- list all data blocks in a format 2 snapshot file -----

def list_format2_blocks(filename):
  if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()
  
  head = snapshot_header(filename)
  format = head.format
  swap = head.swap
  del head
  
  if (format != 2):
    print "not a format 2 snapshot file"
    sys.exit()
            
  print "#   BLOCK   OFFSET   SIZE"
  print "-------------------------"
  
  find_block(filename, format, swap, "XXXX", 0, only_list_blocks=True)
  
  print "-------------------------"
