# Returns a list of bam files from the nominated directory
setwd("/data/home/apattison/ShinyApps/andrew/")
find_bam_files <- function(file_path) {
  if(file.exists(paste0(file_path,'/','plotter-config.json'))){
    json <- fromJSON(paste0(file_path,'/', "plotter-config.json"))
    bam_files <- json$samples
  }
  else{
    bam_files <- list.files(paste (file_path), pattern = '*.bam$')
  }
  return(bam_files)
}

# Returns a list of gff files from the nominated directory
find_gff_files <- function(file_path) {
  if(file.exists(paste0(file_path,'/','plotter-config.json'))){
    json <- fromJSON(paste0(file_path,'/', "plotter-config.json"))
    gff_files <- json$peaks    
  }
  else{
    gff_files <- list.files(paste(file_path), pattern = '*.gff$')    
  }
  return(gff_files)
}

selected_data<- function (data){
  return (paste(data))
}
# Outpus the rows matching the input gene or peak name   
filter_gff_for_rows<- function (gff,names){
  split_names <- strsplit(names, split = " ")
  empty <- data.frame()
  
  for (name in split_names[[1]]){
    index1 <- with(gff, grepl 
                   (ignore.case = T,paste('[=/]{1}',name,'[;/,]',sep=""), gff[,'Information']))
    # Would be nice to find some better regex to get rid of this if statement. 
    # Maybe do this with a GFF parser
    
    index2 <- with(gff, grepl 
                   (ignore.case = T,paste('=',name,'$',sep=""), gff[,'Information']))
    
    output <-gff[index1 | index2, ] 
    if (nrow (output) == 0){
      stop('There are no reads this gene/peak in your selected samples')
    }
    output$input_gene_or_peak <- name
    empty <- rbind(empty, output)
  }
  
  return(empty)
}

# This function gets the poly (A) counts for all given gff rows
get_a_counts <- function(bam_file_path,gff_rows, bam_files, groups, names_from_json){
  reads_report <- data.frame() 
  for (gff_row in 1:nrow(gff_rows)){
    counts_frame <- get_a_counts_gff_row(bam_file_path, gff_rows[gff_row,], 
                                         bam_files, groups, names_from_json)
    if (nrow(counts_frame) == 0){
      next
    }
    counts_frame$gene_or_peak_name <- gff_rows[gff_row, 'input_gene_or_peak']
    
    reads_report <-rbind(reads_report,counts_frame)   
    
  }
  
  return(reads_report)
}

# Parses the BAM files for eah GFF file entry that we are given
get_a_counts_gff_row <- function(bam_file_path,peak, bam_files, groups,names_from_json){
  if (peak[,"Orientation"]== "-"){
    ori <- TRUE    
  }
  else{
    ori <- FALSE
  }
  bam_frame <- data.frame()
  count <- 1
  for (bam_file in bam_files){
    if (substring(bam_file,1,1)=="/"){
      full_file_path <-bam_file
    }
    else{
      full_file_path <-paste(bam_file_path,"/", bam_file, sep ="")
    }
    
    param <- ScanBamParam(what=c('qname','pos','qwidth','strand', 'seq'),
                          tag=c('AN','AD'),flag=scanBamFlag(isMinusStrand=ori) , 
                          which=GRanges(peak [,'Chromosome'],IRanges(
                            peak[,'Peak_Start'], peak[,'Peak_End'] )))
    #Grabs reads overlapping the range specified by the gff row
    result <- scanBam (full_file_path , param = param, isMinusStrand = ori)
    # A check to make sure the adapter bases column is present. 
    #If not, I make a fake one of 0s.
    
    if (length(result [[1]][[6]][[1]])!= length(result [[1]][[5]])){
      result [[1]][[6]][[1]] <- rep(0, length(result [[1]][[5]]))    
    }
    if (length(result [[1]][[6]][[2]])!= length(result [[1]][[5]])){
      result [[1]][[6]][[2]] <- rep(0, length(result [[1]][[5]]))      
    }
    result[[1]][["seq"]] <- as.character(result[[1]][["seq"]])
    if (length(result [[1]][[5]]) == 0){
      stop(paste('There are no reads for at least one peak in ', bam_file))
    }
    
    single_bam_frame <-  data.frame(result) 
    
    colnames(single_bam_frame)<- c("qname", "strand", "pos", 
                                   "width", "sequence", "number_of_as", "number_of_ad_bases")
    #If the read is on the forward strand, add width to pos to obtain 3' end. 
    if (ori == FALSE ){
      single_bam_frame$pos <- single_bam_frame$pos+ single_bam_frame$width
    }
    
    single_bam_frame <- single_bam_frame[single_bam_frame$pos >= 
                                           peak[,'Peak_Start']& 
                                           single_bam_frame$pos <= peak[,'Peak_End'] ,]
    if (nrow(single_bam_frame) == 0){
      next
    }
    if (substring(bam_file,1,1)=="/"){
      single_bam_frame$sample <-  names_from_json$name [names_from_json$bam ==paste(bam_file)]
    }
    else{
      single_bam_frame$sample <- paste(bam_file)
    }
    single_bam_frame$group<- paste("group", groups[count])
    
    bam_frame <- rbind(bam_frame,single_bam_frame)
    count <- count +1
    
  }  
  return(bam_frame)
}

#Strips whitespace out of names
trim <- function (x) gsub("^\\s+|\\s+$", "", x)  

# 
names_string <- function(s_frame, groups, all_reads){
  # The data frame is split by samples here
  to_print <- character()
  for (frame in s_frame){
    #The data frame is split into genes here
    split_peaks <- split(frame ,frame$gene_or_peak_name, drop =T)
    for (peak_frame in split_peaks){
      if (all_reads == T){
        tail_reads <- " "        
      }
      else{
        tail_reads <- " with a poly (A)-tail "
      }      
      if (groups == T){
        str <- paste("The number of reads ",tail_reads,"for ", peak_frame$group[1]," ", 
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
      }
      else{
        str <- paste("The number of reads",tail_reads,"for ",  peak_frame$sample[1]," ",
                     peak_frame$gene_or_peak_name[1], " is: ",nrow(peak_frame),".", "\n", sep ="")
      }
      to_print <- c(to_print, str)
    }
  }
  return(to_print)
}
# Handles overlapping peaks in the gff file 
modify_gff_inplace <- function (gff_file, name) {
  
  saved_gff <-   paste0("/data/home/apattison/ShinyApps/dev/PAT-seq-explorer/gff_dump/",
                        name, ".gff")
  if (file.exists(saved_gff)){
    gff <- read.delim(saved_gff, header=T,
                      comment.char="",stringsAsFactors=F)
    return(gff)
  }
  
  
  start_gff_file <- read.delim(gff_file, header=FALSE,
                               comment.char="",stringsAsFactors=F)
  
  colnames(start_gff_file)<- c('Chromosome', 'Generated_By', 'Feature_Type', 
                               'Peak_Start','Peak_End','-',
                               'Orientation', '--','Information')
  
  plus_frame <- start_gff_file[start_gff_file[,'Orientation'] == '+',]
  plus_frame [,c('Peak_Start', 'Peak_End')] <- plus_frame [,c('Peak_Start', 'Peak_End')]+12
  
  plus_reads <- plus_frame[
    with(plus_frame,order(
      Chromosome,Orientation,Peak_Start)
    ),
    ]
  
  minus_frame <- start_gff_file[start_gff_file[,'Orientation'] == '-',]
  minus_frame [,c('Peak_Start', 'Peak_End')] <- minus_frame [,c('Peak_Start', 'Peak_End')]-12
  
  minus_reads<- minus_frame[
    with(minus_frame,order(
      Chromosome,Peak_Start)
    ),
    ]
  for (row in 1:nrow(plus_reads)){
    if (row == 1){
      next
    }
    if (plus_reads[row, 'Chromosome'] != plus_reads[row-1,'Chromosome']){
      next
    }
    if (plus_reads[row,'Peak_Start'] <= plus_reads[row-1,'Peak_End']){
      plus_reads[row,'Peak_Start'] <- 
        plus_reads[row-1,'Peak_End']+1 
    }
  }
  for (row in 1:nrow(minus_reads)){
    if (row==nrow(minus_reads)){
      next
    }
    if (minus_reads[row, 'Chromosome'] != minus_reads[row+1,'Chromosome']){
      next
      
    }
    if (minus_reads[row,'Peak_End'] >= minus_reads[row+1,'Peak_Start']){
      minus_reads[row,'Peak_End'] <- minus_reads[row+1,'Peak_Start']-1 
    }
  }
  new_frame <- rbind(plus_reads, minus_reads)
  #will new to become getwd
  write.table(x= new_frame, file =saved_gff , append = F,
              quote = F, sep = "\t", row.names = F, col.names = T)
  return(new_frame)
}
# Makes the menas and medians frame shown in info tab of the app
make_means_and_meds_frame <- function (poly_a_counts){
  if (poly_a_counts[1,"group"]== "group NULL"){
    into_samples <- split(poly_a_counts, 
                          list(poly_a_counts$sample, poly_a_counts$gene_or_peak_name))
    mm_frame <- data.frame()
    for (sample in into_samples){
      sample_mean <- mean(sample$number_of_as, na.rm =T)
      sample_median <- median(sample$number_of_as, na.rm =T)
      name <- paste(sample[1, "sample"], sample[1, "gene_or_peak_name"])
      to_bind <- cbind(name,sample_mean, sample_median)
      mm_frame <- rbind(mm_frame,to_bind)
    }
    
  }
  else{
    into_samples <- split(poly_a_counts, poly_a_counts$group)
    mm_frame <- data.frame()
    for (sample in into_samples){
      sample_mean <- mean(sample$number_of_as, na.rm =T)
      sample_median <- median(sample$number_of_as, na.rm =T)
      to_bind <- cbind(sample[1, "group"],sample_mean, sample_median)
      mm_frame <- rbind(mm_frame,to_bind)
    }
  }  
  colnames (mm_frame) <- c("Sample Name", "Mean Poly (A)-Tail Length", "Median Poly (A)-Tail Length")
  return(mm_frame)
}


poly_a_plot <- function (processed_frame, ranges,names, leg = F,group = F){
  new_frame <- processed_frame
  
  if (group == T){
    samples <- split(new_frame, new_frame$group, drop =T)
  }
  else {
    samples <- split(new_frame, new_frame$sample, drop =T)    
  }  
  dummy_ecdf <- ecdf(1:10)
  curve((-1*dummy_ecdf(x)*100)+100, from=ranges[1], to=ranges[2], 
        col="white", xlim=ranges, main= paste(names),
        axes=F, xlab= 'Poly (A) tail length', ylab = 'Percent population (%)', ylim =c(0,100))
  axis(1, pos=0, tick = 25)
  axis(2, pos= 0, at= c(0,25,50,75,100), tick = 25)   
  
  count <- 1  
  
  for (df in samples){
    split_peak <- split(df,df$gene_or_peak_name, drop =T)    
    for(gene_or_peak in split_peak){  
      colours <- rainbow(length(samples)*length(split_peak))
      ecdf_a <- ecdf(gene_or_peak[,"number_of_as"])
      curve((-1*ecdf_a(x)*100)+100, from=ranges[1], to=ranges[2], 
            col=colours[count], xlim=ranges, main= paste(names),
            add=T)
      count <- count +1     
      
    }
    
  }
  # This loop makes a list for the legend. 
  leg_names <- list()
  for (name in names(samples)){
    leg_names <- c(leg_names, paste(name, names(split_peak)))
    
  }
  if (leg ==T){ 
    x_offset <-  length(strsplit(paste(leg_names), "")[[1]])
    legend("topright", 
           legend = leg_names, fill = colours, bty ="n")
  }
}


get_genomic_seq <- function(chr, start, end){  
  seq <- getSeq(Hsapiens,chr,start,end)
  c_seq <-  as.character(seq)
  return(c_seq)
}

igv_plot <- function (processed_frame, ranges,names, leg,group = F, 
                      order_alt = T, alt_cumu_dis,show_poly_a =F, poly_a_pileup=T, gffin){
  
  start <- gffin[1, "Peak_Start"]
  end <- gffin[1,"Peak_End"]+300
  
  #   start <- gffin[1, "Peak_Start"]
  #   end <- gffin[1,"Peak_End"]+400
  
  chr <- gffin[1, "Chromosome"]
  #   in_chr <- as.numeric(as.roman (substring(chr, 4)))
  #   str_chr <- paste0("chr",in_chr)
  #   sequence <- get_genomic_seq(str_chr, start, end)
  #   sequence <- strsplit(sequence , "")
  new_frame <- processed_frame
  
  
  if (gffin[1, "Orientation"] =="-"){
    new_frame <- new_frame[
      with(new_frame,order(
        group, sample, -pos+width, -number_of_as)
      ),
      ]   
  }
  else{
    new_frame <- new_frame[
      with(new_frame,order(
        group, sample, pos+width, -number_of_as)
      ),
      ]
  }
  if (group == T){
    group_status <- "group"
    samples <- split(new_frame, new_frame$group, drop =T)
  }
  else {
    group_status <- "sample"
    samples <- split(new_frame, new_frame$sample, drop =T)    
  }  
  
  count <- list()
  for (sample in samples){
    count <- c(count, 1:nrow(sample))
  }
  count <- as.numeric(unlist(count))
  new_frame$count <- as.numeric(unlist(count))
  
  if (gffin[1, "Orientation"] =="-"){
    new_frame$pos <- new_frame$pos+new_frame$width
    new_frame$bam_read_ends <- new_frame$pos-new_frame$width
    new_frame$poly_a_extension <- new_frame$bam_read_ends -new_frame$number_of_as
  }
  else{
    new_frame$bam_read_ends <- new_frame[,"pos"] - new_frame[,"width"]
    new_frame$poly_a_extension <- new_frame$pos + new_frame$number_of_as
    
  }
  rt <- ggplot(data = new_frame, aes(x= pos, y = count))+
    # scale_x_discrete(labels= sequence[[1]])+
    facet_wrap(as.formula(paste("~", group_status)),ncol = 2)+
    geom_segment(aes(x= pos,xend=bam_read_ends,  y= count ,
                     yend= count, colour = "Alligned reads"))+
    
    xlab(paste(names,"\n","chr", chr,"\n", start,"to", end))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_blank())+
    ylab("Number of reads")+
    theme(axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))+
    scale_colour_manual(values = c("Alligned reads"="green", "Poly (A) tail"="blue" ))
  if (show_poly_a==T){
    if(gffin[1, "Orientation"] =="+"){
      rt = rt+ geom_segment(aes(x= pos, xend=poly_a_extension,  y= count ,
                                yend= count, colour = "Poly (A) tail"))   
      
    }
    else{
      rt = rt+ geom_segment(aes(x=bam_read_ends, xend=poly_a_extension,  y= count ,
                                yend= count, colour = "Poly (A) tail"))        
    }
    rm <- regmatches(gffin[,9], regexpr("id=[^;\\s]+",gffin[,9],perl=T))
    names_list <- gsub(x=rm,pattern="(id=)",
                       replacement="",perl=T)
    rt <- rt + geom_segment(data=gffin, aes_string(x="Peak_Start", xend="Peak_End", y=-1, yend=-1), colour="RED")
    loc <- (gffin$Peak_Start + gffin$Peak_End)/2
        if(length(names_list) > 1){
      ycount <- -1*max(count)/12
      rt <- rt + annotate("text", x = loc, y = ycount, label = names_list)
    }
    
  }
  
  
  return(rt)
  
}



pileup_plot <- function (processed_frame, ranges,names, leg,group = F, 
                         order_alt = T, alt_cumu_dis,show_poly_a =F, poly_a_pileup=T ){
  
  if(order_alt==T){
    
    new_frame <- processed_frame[
      with(processed_frame,order(
        -width, -number_of_as)
      ),
      ]
    ylab <- "Sorted Read Number"
  }
  else{
    new_frame <- processed_frame
    ylab <- "Read Number"
  }
  if (group == T){
    samples <- split(new_frame, new_frame$group, drop =T)
  }
  else {
    samples <- split(new_frame, new_frame$sample, drop =T)    
  }  
  par(bty="l", ps = 10, mar=c(5.1,4.1,4.1,8.1), xpd =T)
  
  if (poly_a_pileup == T ){
    if (length(samples) == 1){
      par(mfrow= c(1,1))
    }
    else if ((length(samples)/2)%%1 == 0){
      par(mfrow= c(as.integer(length(samples)/2),2))
    }
    else{
      par(mfrow= c(as.integer(length(samples)/2)+1,2))        
    }
    for (sample in samples) {
      points <- data.frame(sample$width, sample$number_of_as)
      ymax <- nrow(points)  
      
      count <- 1:ymax
      
      plot(NA,xlim=ranges, ylim = c(0, ymax), xlab= "Number of Bases", ylab = ylab, 
           main= paste(sample[1,'sample']))
      for (i in 1:ymax){
        segments(x0= 0, y0= i,x1= points[i,1], col="purple")
        segments(x0= points[i,1], y0= i,x1= points[i,1] +points[i,2] , col="pink")
        
      }
      
    }
    return()
  }
  ymax <- 0
  for (sample in samples){
    title <- sample[1, 'gene_or_peak_name']
    if (nrow (sample) > ymax){
      ymax <- nrow(sample)
    }
    
  }
  
  if (alt_cumu_dis ==T) {
    dummy_ecdf <- ecdf(1:10)
    curve((-1*dummy_ecdf(x)*100)+100, from=ranges[1], to=ranges[2], 
          col="white", xlim=ranges, main= paste(names),
          axes=F, xlab= "Number of Bases", ylab = 'Percent Population (%)', ylim =c(0,100))
    axis(1, pos=0, tick = 25)
    axis(2, pos= 0, at= c(0,25,50,75,100), tick = 25) 
    count <- 1  
    for (df in samples){
      split_peak <- split(df,df$gene_or_peak_name, drop =T)    
      for(gene_or_peak in split_peak){  
        colours <- rainbow(length(samples)*length(split_peak))
        
        ecdf_a <- ecdf(gene_or_peak[,"width"])
        curve((-1*ecdf_a(x)*100)+100, from=ranges[1], to=ranges[2], 
              col=colours[count], xlim=ranges, main= paste(names),
              add=T)
        count <- count +1     
        
      }
      # This loop makes a list for the legend. 
      leg_names <- list()
      for (name in names(samples)){
        leg_names <- c(leg_names, paste(name, names(split_peak)))
        
      }
      if (leg ==T){ 
        x_offset <-  length(strsplit(paste(leg_names), "")[[1]])
        legend(ranges[2]-30-(x_offset)*2,110 +(length(samples)*-0.8), 
               legend = leg_names, fill = colours, bty ="n")
      }
      
    }
  }
}

gene_expression_plot <- function(processed_bame_frame){
  if (processed_bame_frame[1,"group"] == "group NULL"){
    samples <- split(processed_bame_frame,processed_bame_frame$sample, drop =T)
    xlab <- "Sample"
  }
  else{
    samples <- split(processed_bame_frame,processed_bame_frame$group, drop =T)
    xlab <- "Group"
  }
  df <- data.frame("Sample"= character(), "Count"= numeric())
  for (i in 1:length(samples)){
    
    row <- data.frame(names(samples[i]),nrow(samples[i][[1]]))
    df <- rbind(df,row)
    
  }
  colnames(df) <- c("Sample", "Count")
  gplot <- ggplot(data= df, aes(x=factor(Sample), y = Count))+
    geom_bar(stat =  "identity",colour = "green",, fill ="blue", width =0.5 )+
    xlab(xlab)+
    ylab("Raw number of reads")
  return(gplot)
  
  
}

#help_text (filters out reads that
#were not sequenced completely to the end of the poly (A)-tail
# gets reads that had a 
#  genomic aligment within the selected length
