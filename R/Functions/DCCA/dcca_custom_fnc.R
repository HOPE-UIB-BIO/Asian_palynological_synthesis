##################################################
## R code executing DCCA via external CANOCO.EXE
## from Canoco 4.56 package - just in a specific
## manner for the HOPE project.
## Written by Petr Smilauer, 2019
##################################################

##################################################
# Performs single execution of command-line CANOCO program
#  with specified response and predictor data, using DCCA after
#  square-root transforming. Returns analysis results as
#  a list (described below). On failures calls the stop
#  function, so it might be wise to enclose the call in a try
#  construct to catch such cases
# We need to specify first three arguments to perform the analysis
# The other arguments are optional and adjust behaviour
# Arguments:
#   path   absolute path to the folder where .EXE is located
#          and where also the intermediate results will be stored
#          (so it should be writable)
#   resp   data frame with response data (pollen counts or percentages)
#   pred   data frame with required predictors (all will be used)
#   base   how to name CON, OUT and SOL files (e.g. xx.con, xx.out, xx.sol)
#          Default value "xx"            
#          The 'base' is also used for produced data files: 
#            xx-pred.dta, xx-resp.dta
#   downw  logical - whether to downweight rare species
#          Default value FALSE
#   
# Function returns a list with required information about
#  analysis results. The list contains following named members:
#  eig          numeric vector - eigenvalues for the first four axes
#  tot.inertia  total variation in (transformed) response data
#  turn	        numeric vector with turnover values
#  case.e       numeric matrix with CaseE scores for (up to) first 3 axes
#  case.r       numeric matrix with CaseR scores for first 4 axes
##################################################
execute.can <- function( path, resp, pred, 
                         base="xx", downw=FALSE)
{
  ### (0) Define internally used functions

  ##################################################
  # Saves a data frame into a file in Canoco free format
  #   path is the fully qualified name of the file to create
  #   data is the data frame
  # Function returns FALSE on failure, TRUE otherwise
  ##################################################
  export.can.file <- function( path, data)
  {
    n.rows <- nrow( data)
    n.cols <- ncol( data)
    # Try to create an empty file first
    if(!file.create( path))
      return( FALSE);
    # f.lines will store all output lines
    f.lines <- vector("character")
    f.lines[1] <- "Free format file created from a data frame"
    f.lines[2] <- "FREE"
    f.lines[3] <- paste( n.cols, n.rows)  # space separated
    idx.line <- 4                         # present line index ...
    for( i in 1:n.rows)
    {
      a.line <- ""
      for( j in 1:n.cols)
      {
        a.line <- paste( a.line, data[i,j])
        if((j %% 10) == 0)
        {
          f.lines[idx.line] <- a.line     # store current line state
          a.line <- ""                    # reset to empty
          idx.line <- idx.line + 1        # move to the next line index
        }
      }
      # if we have something in the output buffer ...
      if(nchar(a.line) > 0)
      {
        f.lines[idx.line] <- a.line       # store current line state
        idx.line <- idx.line + 1          # move to the next line index
      }
    }
    # After data, store the names of variables
    a.line <- ""
    col.names <- dimnames( data)[[2]]
    for( i in 1:n.cols)
    {
      a.line <- paste( a.line, 
                     formatC( substring( col.names[i], 1, 8), # up to first 8 chars
                       format="s", width=8, flag="-"), # left-aligned, padded with space
                     sep="")  # no space among the labels
      if((i %% 10) == 0)
      {
        f.lines[idx.line] <- a.line      # store current line
        a.line <- ""                     # reset line to empty
        idx.line <- idx.line + 1         # move to the next line index
      }
    }
    if(nchar( a.line) > 0)               # something left in otput buffer
    { 
      f.lines[idx.line] <- a.line        # store current line
      a.line <- ""                       # set buffer to empty
      idx.line <- idx.line + 1           # progress to next line index
    }
    row.names <- dimnames( data)[[1]]
    for( i in 1:n.rows)
    {
      a.line <- paste( a.line,
                       formatC( substring(row.names[i],1,8),    
                              width=8, format="s", flag="-"), # see above for parameters
                       sep="")
      if((i %% 10) == 0)
      {
        f.lines[idx.line] <- a.line
        a.line <- ""
        idx.line <- idx.line + 1
      }
    }
    if(nchar( a.line) > 0)
    {
      f.lines[idx.line] <- a.line
    }
    # Open the file to be exported for writing
    dta.1 <- file( path, open="wt")
    # Output the lines
    writeLines( f.lines, dta.1)
    # Close the file
    close( dta.1)
    return( TRUE)
  }

  ##################################################
  # Creates CANOCO 4.5x project (.CON) file for a specific
  # type of analysis (DCCA with detrending by segments)
  # pathCON     qualified name of CON file to create
  # nameResp    name of (CANOCO-format) file with response data
  # namePred    name of (CANOCO-format) file with explanatory data
  # base        temporary file name base (for .OUT and .SOL files) 
  #             (default "xx")
  # downwght    logical value: down-weighting of rare species 
  #             (default FALSE)
  # nsegs       number of segments to use for detrending 
  #             (default 26)
  # nresc       number of times to apply non-linear rescaling 
  #             (default 4)
  # rescThresh  threshold for non-linear rescaling 
  #             (default 0.0)
  # Return value: TRUE on success, FALSE otherwise
  ##################################################
  export.CON.file <- function( pathCON, nameResp, namePred, 
                               base="xx", downwght=FALSE,
                               nsegs=26, nresc=4, rescThresh=0.0)
  {
    if(!file.exists( nameResp))
      return( FALSE);
    if(!file.exists( namePred))
      return( FALSE);
    if(!file.create( pathCON))
      return (FALSE);

    f.lines <- vector("character")
    f.lines[ 1]  <- "     2"
    f.lines[ 2]  <- "     1 =  long dialogue?"
    f.lines[ 3]  <- "  0 = changing maximum sizes?"
    f.lines[ 4]  <- paste( "", nameResp)               # insert space at start
    f.lines[ 5]  <- " S"                               # no covariates
    f.lines[ 6]  <- paste( "", namePred)               # insert space at start
    f.lines[ 7]  <- paste( " ", base, ".out", sep="")  # OUT file, e.g. " xx.out"
    f.lines[ 8]  <- paste( " ", base, ".sol", sep="")  # solution file, e.g. " xx.sol"
    f.lines[ 9]  <- " 8   = analysis number (DCCA)"
    f.lines[10]  <- " 1   = detrending by segments"
    f.lines[11]  <- paste( " ", nsegs,"                  = number of segments", sep="")
    f.lines[12]  <- paste( "  ", nresc, " = rescaling of axes", sep="")
    f.lines[13]  <- paste( "   ", formatC( rescThresh, digits=2, format="f"),
                           "= rescaling threshold")
    f.lines[14]  <- "  2 = no. of axes in spec-env biplot"  
    f.lines[15]  <- "  0 = spec and sample diagnostics"
    f.lines[16]  <- "     0 = sample number to be omitted"
    f.lines[17]  <- "  0 = select/delete of environmental variables"
    f.lines[18]  <- "    0     0  = product of environmental variables"
    f.lines[19]  <- "   2 = square root transformation of species data"
    f.lines[20]  <- "     1.00000 = weight for species  ( noweight=1)"
    f.lines[21]  <- "     1.00000 = weight for  samples ( noweight=1)"
    if(downwght)
      f.lines[22]<- "  1 = weighting of species?"
    else
      f.lines[22]<- "  0 = weighting of species?"
    f.lines[23]  <- "  0 = output of correlations?"
    f.lines[24]  <- "  0  2  0  0  0  0  0  2 = output just CaseR and CaseE"

    # open CON file for writing
    con.file <- file( pathCON, open="wt")
    # write formed lines to it
    writeLines( f.lines, con.file)
    # close it
    close( con.file)  

    return( TRUE)
  }  

  ##################################################
  # Takes an open connection to a text file and reads from it individual
  # lines until it arrives to a line starting with specified text
  # Found line is returned 
  # Passed connection (connFile) must be already open for reading
  # The function returns empty string to indicate a failure
  ##################################################
  read.to.line <- function( connFile, soughtText)
  {
    aLine <- vector("character")
    while( TRUE)
    {
      aLine <- readLines( connFile, n=1)
      if(length( aLine) == 0)
        return( "");            # return empty string on failure
      if(startsWith( aLine, soughtText))
        break;
    }
    return( aLine[1]);
  }

  ##################################################
  # Reads specified number of lines from the passed 
  # connection. Last read line is returned
  ##################################################
  read.n.lines <- function( connFile, numLines)
  {
    aLine <- readLines( connFile, n=numLines)
    nL <- length( aLine)
    if(nL == 0)
      return( "")
    return( aLine[nL]);
  }

  ##################################################
  # Takes a character string and parses it with TAB separator, 
  #  returning specified range of string entries
  ##################################################
  parse.line.by.tabs <- function( aLine, idxFrom, idxTo)
  {
    xy <- strsplit( aLine, "\t", fixed=T) [[1]]
    xy[idxFrom:idxTo]
  }


  ######################################################
  ### (1) Create data files
  path.len <- nchar( path)
  if(substr( path, path.len, path.len) != "/")  # append slash to path, if not present
   path <- paste( path, "/", sep="")
  # fully qualified name of response data file   
  resp.fname <- paste( path, base, "-resp.dta", sep="")
  # create response data file
  if(export.can.file( resp.fname, resp) == FALSE)
    stop("Cannot export response data into ", resp.fname, call.=FALSE)
  
  # fully qualified name of  predictor data file
  pred.fname <- paste( path, base, "-pred.dta", sep="")
  # create predictor data file
  if(export.can.file( pred.fname, pred) == FALSE)
  {
    file.remove( resp.fname)
    stop( "Cannot export predictor data into ", pred.fname, call.=FALSE)
  }

  ### (2) Create CON file
  con.fname <- paste( path, base, ".con", sep="")
  if(export.CON.file( con.fname, resp.fname, pred.fname, base, downw)==FALSE)
  {
    file.remove( resp.fname);
    file.remove( pred.fname);
    stop( "Cannot produce CON file into ", con.fname, call.=FALSE)
  }

  ### (3) Execute CANOCO command line
  WD <- getwd()
  setwd( path)
  system2( "Canoco.exe", args=con.fname)

  ### (4) Parse OUT file
  out.fname <- paste( base, ".out", sep="")
  if(file.exists( out.fname) == FALSE)
  {
    file.remove( resp.fname);
    file.remove( pred.fname);
    file.remove( con.fname);
    stop("Cannot find produced file ", out.fname, call.=FALSE)
  }
  out.file <- file( out.fname, open="r")
  # first determine the number of active sample:
  aLine <- read.to.line( out.file, " No. of active  samples:");
  nrows <- as.numeric( (strsplit( aLine, ":", fixed=T) [[1]])[2])
  nrows.tot <- dim( resp)[1]
  if((nrows < 1) | (nrows > nrows.tot))
  {
    file.remove( resp.fname);
    file.remove( pred.fname);
    file.remove( con.fname);
    stop("Wrong count of active samples ", nrows, " in output", call.=FALSE)
  }
  aLine <- read.to.line( out.file, " Eigenvalues  ");
  if(nchar(aLine) < 20)
  {
    file.remove( resp.fname, pred.fname, con.fname);
    # file.remove( out.fname);
    stop("Cannot parse produced file ", out.fname, call.=FALSE)
  }
  # Retrieve eigenvalues and total inertia
  eigs <- as.numeric( parse.line.by.tabs( aLine, 2, 6))
  
  tot.inertia <- eigs[5]
  eigs <- eigs[1:4]
  
  # Retrieve and parse turnover size
  aLine <- read.n.lines( out.file, 1)
  turns <- as.numeric( parse.line.by.tabs( aLine, 2, 5))
  # close OUT file
  close( out.file)
   
  
  ### (5) Parse SOL file
  sol.fname <- paste( base, ".sol", sep="")
  if(file.exists( sol.fname) == FALSE)
  {
    file.remove( resp.fname, pred.fname, con.fname, out.fname);
    stop("Cannot find produced file ", sol.fname, call.=FALSE)
  }
  sol.file <- file( sol.fname, open="r")
  aLine <- read.to.line( sol.file, " Samp: Sample scores");
  if(nchar(aLine) < 20)
  {
    file.remove( resp.fname, pred.fname, con.fname, out.fname);
    stop("Cannot parse produced file [1] ", sol.fname, call.=FALSE)
  }
  
  # read through another five lines
  aLine <- read.n.lines( sol.file, 5) 
  # prepare storage full of NA values
  case.r <- matrix( as.vector( rep(NA, nrows.tot*4), mode="numeric"), 
                    ncol=4)
  for( idx in 1:nrows)
  { # read single line
    aLine <- read.n.lines( sol.file, 1)
    # separate its fields
    str.entries <- strsplit( aLine, "\t", fixed=T)[[1]]
    # identify case index
    idxCase <- as.numeric( str.entries[1])
    if((idxCase < 1) | (idxCase > nrows.tot) | (length(str.entries) < 8))
    {
      file.remove( resp.fname, pred.fname, con.fname, out.fname);
      stop("Error parsing CaseR scores for entry ", idx, call.=FALSE)
    }
    # store four CaseR scores
    case.r[idxCase,] <- as.numeric( str.entries[3:6]) 
  }
  # Assign row and column names
  dimnames( case.r) <- list( dimnames( resp)[[1]], 
                             c("Axis 1","Axis 2","Axis 3","Axis 4"))
  
  # Now read the second set of case scores
  aLine <- read.to.line( sol.file, " SamE: Sample scores")
  if(nchar(aLine) < 20)
  {
    file.remove( resp.fname, pred.fname, con.fname, out.fname);
    stop("Cannot parse [2] produced file ", sol.fname, call.=FALSE)
  }
  aLine <- read.n.lines( sol.file, 5)
  # prepare storage full of NA values
  case.e <- matrix( as.vector( rep( NA, nrows.tot*4), mode="numeric"), 
                    ncol=4)
  for( idx in 1:nrows)  
  {
    # read one line
    aLine <- read.n.lines( sol.file, 1)
    # separate its fields
    str.entries <- strsplit( aLine, "\t", fixed=T)[[1]]
    # identify case index
    idxCase <- as.numeric( str.entries[1])
    if((idxCase < 1) | (idxCase > nrows.tot) | (length(str.entries) < 7))
    {
      file.remove( resp.fname, pred.fname, con.fname, out.fname);
      stop("Error parsing CaseE scores for entry ", idx, call.=FALSE)
    }
    # store four CaseE scores
    case.e[idxCase,] <- as.numeric( str.entries[3:6]) 
  }
  # Assing row and column names
  dimnames( case.e) <- list( dimnames( resp)[[1]], 
                             c("Axis 1","Axis 2","Axis 3","Axis 4"))
  # close SOL file
  close( sol.file)


  setwd( WD)          # undo change of working directory of R

  ### (6) Return results
  list( eig = eigs, tot.inertia=tot.inertia, turn = turns, 
        case.e = case.e, case.r = case.r)
}


##################################################
# This is a driver function that finds a specific DCCA model
# depending on the comparison of first constrained and first
# unconstrained eigenvalue
# The single numerical predictor (depth or age) is tried first
#  as linear, then as poly(x,2), and finally as poly(x,3)
#  until L(1) > L(unc1)
# You need to specify first three arguments to perform the analysis
# The other two arguments are optional and adjust behaviour
# Function returns a list containing required information about
#  chosen DCCA type and its results
# Arguments:
#   path      is absolute path to the folder where CANOCO.EXE is located
#             and where also the intermediate results will be stored
#             (so it should be writable)
#   resp      data frame with response data (pollen counts or percentages)
#   pred.var  numeric vector with the predictor values (depth or age)
#             Its length must be identical with the count of rows in resp 
#             data frame
#   base      how to name the CON, OUT and SOL files (e.g. xx.con, xx.out, 
#             xx.sol)
#             Also used for produced data files: xx-pred.dta, xx-resp.dta
#   downw     logical, whether to downweight rare species or not
#   
# Return value is a list with following components:
#  eig          numeric vector with the eigenvalues for first four axes
#  tot.inertia  total variation in response data
#  turn	        numeric vector with turnover values for first four axes
#  case.e       numeric matrix with CaseE scores
#  case.r       numeric matrix with CaseR scores
#  degree       indicates how the pred.var was used: 1 - as a linear predictor 
#               2 or 3 - as a corresponding polynomial
##################################################
select.DCCA.can <- function( path, resp, pred.var, 
                             base="xx", downw=FALSE)
{
  # Create data frame with a linear form of predictor
  pred.df <- data.frame( x=pred.var)
  dimnames( pred.df) <- list( dimnames( resp)[[1]], c("Predictr"))
  # Execute CANOCO
  res.1 <- execute.can( path, resp, pred.df, base, downw);
  if(length( res.1$eig) < 4)
    stop("Error executing DCCA with a linear predictor effect")
  if(res.1$eig[1] > res.1$eig[2])  # appropriate result
  {
    res.1$degree <- 1
    return( res.1)   
  }
  
  # Failed, so create 2nd-order polynomial
  #  Make a copy of predictor ...
  pred.df[,2] <- pred.df[,1]  
  # ... and replace both columns with an orthogonal polynomial
  pred.df[,1:2] <- poly( pred.df[,1], 2)
  names( pred.df) <- c("Poly.1", "Poly.2")
  # ... re-execute DCCA
  res.1 <- execute.can( path, resp, pred.df, base, downw)
  if(length( res.1$eig) < 4)
    stop("Error executing DCCA with a poly(x,2) predictor effect")
  if(res.1$eig[1] > res.1$eig[3])  # success
  {
     res.1$degree <- 2
     return( res.1)
  }
  
  # No help, so create 3rd-order polynomial
  pred.df[,3] <- pred.var
  pred.df[,1:3] <- poly( pred.var, 3)
  names( pred.df) <- c("Poly.1", "Poly.2", "Poly.3")
  # ... re-execute DCCA
  res.1 <- execute.can( path, resp, pred.df, base, downw)
  if(length( res.1$eig) < 4)
    stop("Error executing DCCA with a poly(x,3) predictor effect")
  if(res.1$eig[1] < res.1$eig[4])
    warning("Failed to get eigenvalue(Ax1) greater than the 1st unconstrained eigenvalue")
  res.1$degree <- 3
  res.1
}
