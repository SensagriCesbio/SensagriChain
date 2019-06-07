import csv

# CesBIO 2018 #
# Allow to get date from date file

def dateconversion(date):
  monthlength = [31,29,31,30,31,30,31,31,30,31,30,31]; # Attention 2016 is a bisextile year !!!
  monthlabel  = ["Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]  
  d = int(date)%100 
  m = (int(date)/100)%100
  y = (int(date)/10000)%100
  #print date,d,m,y
  label = str(d) + "-" + monthlabel[(m - 1)] 
  mt = 0
  for i in range(m-1):
    mt = mt + monthlength[i]
  stamp = (y-16)*366 + mt + (d-1)
  return [stamp, label]

def stamp2string(d):
  monthlength = [31,29,31,30,31,30,31,31,30,31,30,31]; # Attention 2016 is a bisextile year !!!
  monthlabel  = ["Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]  
  d = int(d)
  if d < 0:
      dt = d
      i = 11
      while(dt < 0):
          day = dt 
          dt = dt + monthlength[i%12]  
          i -= 1
      stamp = str(monthlength[(i-1)%12] + day) + "-" +monthlabel[(i+1)%12]
  elif d == 0:
      stamp = "1-Jan"
  else:
      dt = d
      i = 0
      while(dt > 0):
          day = dt 
          dt = dt - monthlength[i%12]  
          i += 1
      stamp = str(day) + "-" +monthlabel[(i-1)%12]
      
  return stamp

def GetDates(datefile):
  with open(datefile) as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)

  DoY = [i[0] for i in d]
  DoYlabel = []
  date = []
  if int(DoY[0]) > 20000000:  
      for d in DoY:
          dc = dateconversion(d)
          date.append(dc[0])
          DoYlabel.append(dc[1])
  else:
      for d in DoY:
          date.append(int(d))
          DoYlabel.append(stamp2string(d))

  return date,DoYlabel

def GetDatesSentinel(datefile):
  with open(datefile) as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)

  DoY = [i[0][17:25] for i in d]
  DoYlabel = []
  date = []

  if int(DoY[0]) > 20000000:  
      for d in DoY:
          dc = dateconversion(d)
          date.append(dc[0])
          DoYlabel.append(dc[1])
  else:
      for d in DoY:
          date.append(int(d))
          DoYlabel.append(stamp2string(d))

  return date,DoYlabel


