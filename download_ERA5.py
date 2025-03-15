# Make changes in the dates only
year = "2016"
month = "02"
days = ["21","22","23"]
###
# Import cdsapi and create a Client instance
import cdsapi
c = cdsapi.Client()
c.retrieve("reanalysis-era5-pressure-levels", {
        "product_type":   "reanalysis",
        "area":           "60.00/-20.00/40.00/20.00",
        "variable":       ["z","t","r"],
        "pressure_level": ["1","2","3","5","7","10","20","30","50","70","100",
                           "125","150","175","200","225","250","300","350","400",
                           "450","500","550","600","650","700","750","775","800",
                           "825","850","875","900","925","950","975","1000"],
        "year":           year,
        "month":          month,
        "day":            days,
        "time":           ["00","01","02","03","04","05","06","07","08","09","10","11",
                           "12","13","14","15","16","17","18","19","20","21","22","23"],
    }, "PRES_SC_"+year+"_"+month+"_"+days[0]+"-"+days[-1]+".grb")


# Import cdsapi and create a Client instance
import cdsapi
c = cdsapi.Client()
c.retrieve("reanalysis-era5-pressure-levels", {
        "product_type":   "reanalysis",
        "area":           "60.00/-20.00/40.00/20.00",
        "variable":       ["u","v"],
        "pressure_level": ["1","2","3","5","7","10","20","30","50","70","100",
                           "125","150","175","200","225","250","300","350","400",
                           "450","500","550","600","650","700","750","775","800",
                           "825","850","875","900","925","950","975","1000"],
        "year":           year,
        "month":          month,
        "day":            days,
        "time":           ["00","01","02","03","04","05","06","07","08","09","10","11",
                           "12","13","14","15","16","17","18","19","20","21","22","23"],
    }, "PRES_UVW_"+year+"_"+month+"_"+days[0]+"-"+days[-1]+".grb")

# Import cdsapi and create a Client instance
import cdsapi
c = cdsapi.Client()
c.retrieve("reanalysis-era5-single-levels", {
        "product_type":   "reanalysis",
        "area":           "60.00/-20.00/40.00/20.00",
        "variable":       ["10u","10v","2t","2d","msl","sp","sst","skt",
                           "stl1","stl2","stl3","stl4","slt","swvl1","swvl2","swvl3","swvl4",
                           "sd","rsn","lsm","ci"],
        "year":           year,
        "month":          month,
        "day":            days,
        "time":           ["00","01","02","03","04","05","06","07","08","09","10","11",
                           "12","13","14","15","16","17","18","19","20","21","22","23"],
    }, "SFC_"+year+"_"+month+"_"+days[0]+"-"+days[-1]+".grb")

