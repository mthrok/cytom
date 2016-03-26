#
# ImmPort FCS conversion program: FCSTrans
# version: 1.3 - accuri data support added
# Authors: Yue Liu and Yu "Max" Qian
# Contact: yu.qian@utsouthwestern.edu or yliu0@yahoo.com
# Usage: 1)	Load the flowCore package and other relevant packages in R if necessary
#	 2)	source(FCSTrans.R)
#	 3)	ipconvert(test.fcs)
#	        or ipconvert(./FCSfoldername)
#
# FCSTrans has been published at: http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22037/abstract
#
# FCSTrans can be cited as:
# Qian, Y., Liu, Y., Campbell, J., Thomson, E., Kong, Y. M. and Scheuermann, R. H. (2012), FCSTrans: An open source software system for FCS file conversion and data transformation. Cytometry. doi: 10.1002/cyto.a.22037

library(marray)
library(flowCore)

# set output to 0 when input is less than cutoff value
ipfloor <- function (x, cutoff = 0, target = 0) {
    y = x
    if (x <= cutoff) y = target
    y
}

# set output to 0 when input is less than cutoff value
ipceil <- function (x, cutoff = 0, target = 0) {
    y = x
    if (x >= cutoff) y = target
    y
}

# calculation core of iplogicle
iplogicore <- function (x, w, r, d, scale) {
    tol = .Machine$double.eps^0.8
    maxit = as.integer(5000)
    d = d * log(10)
    scale = scale / d
    p = if (w == 0) 
        1
    else uniroot(function(p) -w + 2 * p * log(p)/(p + 1), c(.Machine$double.eps, 
        2 * (w + d)))$root
    a = r * exp(-(d - w))
    b = 1
    c = r * exp(-(d - w)) * p^2
    d = 1/p
    f = a * (p^2 - 1)
    y = .Call("biexponential_transform", as.double(x), a, b, c, d, f, w, tol, maxit)
    y = sapply(y * scale, ipfloor)
    y
}

# function for calculating w 
iplogiclew <- function (w, cutoff = -111, r = 262144, d = 4.5, scale = 1) {
    if (w > d) 
        w = d
    y = iplogicore(cutoff, w, r, d, scale) - .Machine$double.eps^0.6
    y
}

# immport logicle function - convert fluorescent marker values to channel output
iplogicle <- function (x, r=262144, d=4.5, range=4096, cutoff=-111, w=-1) {
    if (w > d) 
        stop("Negative range decades must be smaller than total number of decades")
    if (w < 0)
        w = uniroot(iplogiclew, c(0, d), cutoff=cutoff)$root
    print(paste("params: r=", r, "d=", d, "range=", range, "cutoff=", cutoff, "w=", w))
    y = iplogicore(x, w, r, d, range)
    y
}

# Convert fluorescent values to channel output using log transformation
# logicle transformation after compensation
iplog <- function(x) {
    x = sapply(x, ipfloor, cutoff=1, target=1)
    y = 1024 * log10(x) - 488.6
    y
}

# immport linear function - convert scatter values to channel output
# linear transformation
ipscatter <- function (x, channelrange=262144) {
    y = 4095.0 * x / channelrange
    y = sapply(y, ipfloor)
    y = sapply(y, ipceil, cutoff=4095, target=4095)
    y
}

# immport time function - convert time values to channel output
# linear transformation
iptime <- function (x, channelrange) {
    # use simple cutoff for now
    y = sapply(x, ipfloor)
    y
}

# get marker type
getMarkerType <- function(name) {
    type = ""
    prefix2 = toupper(substr(name, 1, 2))
    prefix3 = toupper(substr(name, 1, 3))
    prefix4 = toupper(substr(name, 1, 4))
    if (prefix2 == "FS" || prefix2 == "SS") {
        type = "SCATTER"
    } else if (prefix3 == "FSC" || prefix3 == "SSC") {
        type = "SCATTER"
    } else if (prefix4 == "TIME") {
        type = "TIME"
    } else {
        pieces = unlist(strsplit(name, "-"))
        if (toupper(pieces[length(pieces)]) == "A") {
            type = "FLUO_AREA"
        } else {
            type = "FLUO_NON_AREA"
        }
    }
    print(paste("marker:", name, ", type:", type))
    type
}

scaleData <- function(data, channelrange=0) {
    datamax = range(data)[2]  # range() returns [min, max]
    if (datamax > channelrange) {
        channelrange = datamax
    }
    #if (channelrange == 0) {
    #    channelrange = range(data)[2]
    #}
    data = 262144 * data / channelrange
    data
}

isAccuriData <- function(keydata) {
    isTRUE(as.character(keydata$"$CYT") == "Accuri C6")
}

# immport convert function for Accuri FCS data file
# @param fcs    flowFrame object
convertAccuriFcs <- function(fcs) {
    debug = TRUE
    keywords = keyword(fcs)
    fcs.exprs = exprs(fcs)
    
    #apply compensation if available
    #comment the following line if no compensation is necessary
    spill = keywords$'$SPILLOVER'  

    if (is.null(spill) == FALSE) {
        chunks = unlist(strsplit(spill, ','))
        # first chunks is spill matrix size
        msize = as.integer(chunks[1])
        # the next msize of chunks are fluorescent channel numbers
        mchannels = as.integer(chunks[2:(msize+1)])
        # the rest of chunks are the actual spill matrix
        mdata = as.double(chunks[(msize+2):length(chunks)])
        spillmatrix = matrix(mdata, nrow=msize, byrow=T)
        # apply compensation matrix (inverse of spillmatrix) to data
        compexp = fcs.exprs[, mchannels] %*% solve(spillmatrix)
        # put compensated data back in 
        for (i in 1:length(mchannels)) { fcs.exprs[,mchannels[i]] = compexp[,i] }
    }
    
    #process fcs expression data by marker type
    fcs.channel = NULL
    markers = colnames(fcs)
    
    if (debug) print("loop through markers")
    for (i in 1:length(markers)){
        markertype = getMarkerType(markers[i])
        channelrange = 16777215
        channeldecade = 7.224719870049579
        if (markertype == "SCATTER") {
            # 2x zoom on scatter channels for Accuri data
            channel = ipscatter(fcs.exprs[, i], channelrange/1.0)
        } else if (markertype == "TIME") {
            channel = iptime(fcs.exprs[, i])
        } else {
            # apply logicle transformation on fluorescent channels
            channel = iplogicle(fcs.exprs[, i], channelrange, channeldecade)
        }
        fcs.channel = cbind(fcs.channel, round(channel))
    }
    colnames(fcs.channel) = markers
    fcs.channel
}

# immport convert function - convert flow cytometry values to channel output
# iterate columns name and treat data in three categories:
#   scatter      linear 
#   time         linear
#   fluorescent  logicle
# example:
#   > ad008 = read.FCS('ad008.fcs', transformation=F)
#   > ad008_ipc = convertfcs(ad008r)
# @param fcs    flowFrame object
convertfcs <- function(fcs) {
    debug = TRUE
    #check file type and FCS version
    if (class(fcs)[1] != "flowFrame") {
        print("convertfcs requires flowFrame object as input")
        return(FALSE)
    }
    keywords = keyword(fcs)
    if (debug) print(paste("FCS version:", keywords$FCSversion))
    if (debug) print(paste("File data type:", keywords['$DATATYPE']))
    if (keywords$FCSversion == "2" || keywords$FCSversion == "3" ) {
        datatype = unlist(keywords['$DATATYPE']) #check data type
        if (datatype == 'F') {
            #apply compensation if available
            spill = keyword(fcs)$SPILL
            
            if (debug) print("check spill")
            if (is.null(spill) == FALSE) {
                tryCatch({ fcs = compensate(fcs, spill) }, error = function(ex) { str(ex); })
            }
            
            if (debug) print("get expression")
            #process fcs expression data by marker type
            fcs.exprs = exprs(fcs)
            fcs.channel = NULL
            markers = colnames(fcs)
            
            if (debug) print("loop through markers")
            for (i in 1:length(markers)){
                markertype = getMarkerType(markers[i])
                rangekeyword = paste("$P", i, "R", sep="")
                #if (debug) print(paste("  range keyword:", rangekeyword))
                if (debug) print(paste("  range value:", keywords[rangekeyword]))
                channelrange = as.numeric(keywords[rangekeyword])
                if (markertype == "SCATTER") {
                    channel = ipscatter(scaleData(fcs.exprs[, i], channelrange))
                } else if (markertype == "TIME") {
                    channel = iptime(fcs.exprs[, i])
                } else {
                    # apply logicle transformation on fluorescent channels
                    channel = iplogicle(scaleData(fcs.exprs[, i], channelrange))
                }
                fcs.channel = cbind(fcs.channel, round(channel))
            }
            colnames(fcs.channel) = markers
        } else if (datatype == 'I') {
            fcs.channel = exprs(fcs)
        } else {
            print(paste("Data type", datatype, "in FCS 3 is not supported"))
            fcs.channel = FALSE
        }
    } else {
        print(paste("FCS version", keyword(fcs)$FCSversion, "is not supported"))
        fcs.channel = FALSE
    }
    fcs.channel
}

# convert fcs into channel output and save result
# @param fcsfile
# @param verbose       Prints file and folder names when TRUE
# @param overwrite     Convert and overwrite existing channel files
# @param renameSource  Rename fcs file to have .fcs suffix
# @param convert       Flag for convert FCS file
# @param savekeys      Flag for save keywords to file
#
# example:
#   > convertfcs('~/data/ad008')
#   This will convert ad008.fcs and create ad008_ip_channel.txt
processfcsfile <- function(fcsfile, verbose=F, overwrite=F,
    renameSource=F, convert=T, savekeys=T) {
    
    pieces = unlist(strsplit(fcsfile, .Platform$file.sep))
    filename = pieces[length(pieces)]
    filepieces = unlist(strsplit(filename, '\\.'))
    
    #replace .fcs with .txt; append .txt if not ending in .fcs
    if (filepieces[length(filepieces)] == 'fcs') {
        filepieces[length(filepieces)] = 'txt'
    } else {
        filepieces[length(filepieces)+1] = 'txt'
    }
    
    pieces[length(pieces)] = paste(filepieces, collapse = '.')
    channelfile = paste(pieces, collapse = .Platform$file.sep)
    pieces[length(pieces)] = paste(c(filename, 'keylist.txt'), collapse='.')
    keylistfile = paste(pieces, collapse = .Platform$file.sep)
    #print(paste("    keylist file:", keylistfile))
    
    doconvert = T
    if (file.exists(channelfile)) {
        doconvert = overwrite
        if (!doconvert) {
            if (verbose) print(paste("      Skipping", filename))
        }
    }
    if (doconvert) {
        fcs <- read.FCS(fcsfile, transformation=F)
        pdata <- keyword(fcs)  #extract keywords from fcs file
        if (convert) {
            if (verbose) print(paste("    Converting", fcsfile))
            if (isAccuriData(pdata)) {
                cdata <- convertAccuriFcs(fcs)
            } else {
                cdata <- convertfcs(fcs)
            }
            write.table(cdata, file=channelfile, quote=F, row.names=F,
                col.names=T, sep='\t', append=F)
            if (verbose) print(paste("      Converted to", channelfile))
        }
        if (savekeys) {
            write.table(as.matrix(pdata), file=keylistfile, quote=F, row.names=T,
                col.names=F, sep='=', append=F)
            if (verbose) print(paste("      Keywords saved to", keylistfile))        
        }
    }
    
    if(renameSource) {
        filepieces = unlist(strsplit(filename, '\\.'))
        if (filepieces[length(filepieces)] != 'fcs') {
            newfcsfile = paste(fcsfile, 'fcs', sep='.')
            print(paste('    ', fcsfile, 'renamed to', newfcsfile))
            file.rename(fcsfile, newfcsfile)
        }
    }
}

# convert fcs files
# @param entry         Can be fcs object, fcs file or folder name
# @param verbose       Prints file and folder names when TRUE
# @param overwrite     Convert and overwrite existing channel files
# @param renameSource  Rename fcs file to have .fcs suffix
#
# exmaple:
#   > ipconvert('~/data/test')
#   > ipconvert(c('/tmp/a.fcs', 't101.fcs', '/data')
ipconvert <- function(entry, verbose=F, overwrite=F, renameSource=F, level=1) {
    entryclass = class(entry)
    if (is.vector(entry) && length(entry) > 1) {
        if (verbose) print (paste(level, ":", "vector - size", length(entry)))
        sapply(entry, ipconvert, verbose=verbose, overwrite=overwrite,
            renameSource=renameSource, level=level+1)
    } else if (entryclass == 'flowFrame') {
        if (verbose) print (paste(level, ":", "flowFrame"))
        #convertfcs(entry)
    } else if (entryclass == 'character') {
        if (file.exists(entry)) {
            if (file.info(entry)$isdir) {
                if (verbose) print (paste(level, ":", "Folder :", entry))
                #get folder name without ending separator
                pattern = paste(.Platform$file.sep, '$', sep='')
                entry = sub(pattern, '', entry)
                files = list.files(entry)
                if (length(files) > 0) {
                    files = paste(entry, files, sep=.Platform$file.sep)
                    sapply(files, ipconvert, level=level+1, verbose=verbose,
                        overwrite=overwrite, renameSource=renameSource)
                }
            } else {
                isvalid = F
                tryCatch({
                    isvalid = isFCSfile(entry)
                }, error = function(ex) {
                    if (verbose) print (paste("    ! Error in isFCSfile", ex))
                })
                if (isvalid) {
                    if (verbose) print (paste(level, ":", "FCS file :", entry))
                    processfcsfile(entry, verbose=verbose, overwrite=overwrite,
                        renameSource=renameSource)
                } else {
                    if (verbose) print (paste(level, ":", "Not FCS file :", entry))
                }
            }
        } else {
            if (verbose) print(paste(level, ":", entry, "does not exist"))
        }
    } else {
        if (verbose) print(paste(level, ":", entry, "is not supported"))
    }
    
    result = ''
    if (level == 1) {
        result = '--- ipconvert summary goes here ---'
    }
}
