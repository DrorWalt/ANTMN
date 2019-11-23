# Author: Dror Walter
# Creating Topic Model Networks (ANTMN Method): Supplementary code
# To Cite:
# Walter D. & Ophir Y. (2019) News Frame Analysis: an Inductive Mixed Method Computational Approach. Communication Methods and Measures. http://dx.doi.org/10.1080/19312458.2019.1639145
# Online published (July 2019)
################################


################## Main Function
## LDAobject is the output of the topicmodels package LDA command
## deleted_topics can be left empty or can be passed a vector of topic numbers to be ommited from the network.
## topic_names can be left empty or passed a vector of strings, each a topic name.
## To calculate topic size data from the pre-process termdoc matrix needs to be implements. 
#### For more elaboration on this option see sample code. As a default topic size is ignored
## save_filename can be left empty or passed a string for the filename to save. Notice that ".graphml" will be added automatically to filename
## *.graphml, while not associated automatically with Gephi, can be opened in Gephi using the open graph menu.

network_from_LDA<-function(LDAobject,deleted_topics=c(),topic_names=c(),save_filename="",topic_size=c()) {
  # Importing needed packages
  require(lsa) # for cosine similarity calculation
  require(dplyr) # general utility
  require(igraph) # for graph/network managment and output

  print("Importing model")
  
  # first extract the theta matrix form the topicmodel object
  theta<-LDAobject@gamma
  # adding names for culumns based on k
  colnames(theta)<-c(1:LDAobject@k)
  
  # claculate the adjacency matrix using cosine similarity on the theta matrix
  mycosine<-cosine(as.matrix(theta))
  colnames(mycosine)<-colnames(theta)
  rownames(mycosine)<-colnames(theta)
  
  # Convert to network - undirected, weighted, no diagonal
  
  print("Creating graph")
  
  topmodnet<-graph.adjacency(mycosine,mode="undirected",weighted=T,diag=F,add.colnames="label") # Assign colnames
  # add topicnames as name attribute of node - importend from prepare meta data in previous lines
  if (length(topic_names)>0) {
    print("Topic names added")
      V(topmodnet)$name<-topic_names
      } 
  # add sizes if passed to funciton
  if (length(topic_size)>0) {
    print("Topic sizes added")
    V(topmodnet)$topic_size<-topic_size
  }
  newg<-topmodnet
  
  # delete 'garbage' topics
  if (length(deleted_topics)>0) {
    print("Deleting requested topics")
    
  newg<-delete_vertices(topmodnet, deleted_topics)
  }
  
  # run community detection and attach as node attribute
  print("Calculating communities")
  
  mylouvain<-(cluster_louvain(newg)) 
  mywalktrap<-(cluster_walktrap(newg)) 
  myspinglass<-(cluster_spinglass(newg)) 
  myfastgreed<-(cluster_fast_greedy(newg)) 
  myeigen<-(cluster_leading_eigen(newg)) 

  V(newg)$louvain<-mylouvain$membership 
  V(newg)$walktrap<-mywalktrap$membership 
  V(newg)$spinglass<-myspinglass$membership 
  V(newg)$fastgreed<-myfastgreed$membership 
  V(newg)$eigen<-myeigen$membership 
  
  # if filename is passsed - saving object to graphml object. Can be opened with Gephi.
  if (nchar(save_filename)>0) {
    print("Writing graph")
    write.graph(newg,paste0(save_filename,".graphml"),format="graphml")
  }
  
  # graph is returned as object
  return(newg)
}

#Example:
mynewnet<-network_from_LDA(LDAobject=LDA.66,
                           deleted_topics=c(5,6,11,12,20,27,37),
                           save_filename="my_net_file")

################## Sample Code
library(topicmodels) #used for topic model estimation
library(ldatuning) # used for K selection
library(quanteda) # for text handling and pre-processing
library(dplyr) # general utility
library(xlsx) #writing excel files
library(parallel) # used for parallel computing when running models

# Data Prep
## "data" is a dataframe object in which each row is a document, and column cosist of the document data.
## The column containing the text for analysis NUST be titled "text".
## we are adding index to be able to match documents across data types. 
# In addition, if the data dataframe has only on column it forces the data to be a dataframe when manipulating it.
data$index<-seq(1:nrow(data))

### removing extremely short documents.
removed_short<-subset(data,nchar(as.character(data$text))<600)
data2<-subset(data,!nchar(as.character(data$text))<600)

### removing duplicate documents
removed_df<-data2[duplicated(data2$text),]
data3 <- data2[!duplicated(data2$text),]

### Text pre-processing
##### import data to quanteda format
mycorpus <- corpus(data3)
##### using quanteda stopwords, with single letters as well
stopwords_and_single<-c(stopwords("english"),LETTERS,letters)
##### preparing dfm obeject. No stemming due to its impact on topic quality
dfm_counts <- dfm(mycorpus,tolower = TRUE, remove_punct = TRUE,remove_numbers=TRUE, 
                  remove = stopwords_and_single,stem = FALSE,
                  remove_separators=TRUE) 
##### trimming tokens too common or too rare to imporve efficiency of modeling
dfm_counts2<-dfm_trim(dfm_counts, max_docfreq = 0.99, min_docfreq=0.005,docfreq_type="prop")
##### converting to LDA ready object
dtm_lda <- convert(dfm_counts2, to = "topicmodels")

## Selecting the appropriate number of topics
result <- FindTopicsNumber(
  dtm_lda,
  topics = c(1:10 * 10, 1:4 * 20 + 100, 0:2 * 50 + 200, seq(from = 80, to = 100, by = 1)),
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  mc.cores = detectCores(),
  verbose = TRUE
)

FindTopicsNumber_plot(result) # Based on the plot, 66 seems the most efficent model. 

# running the model
LDA.66<- LDA(dtm_lda, k=66, method = "Gibbs")

# extracting excel matrices for topic interpretation
LDAfit<-LDA.66

mybeta<-data.frame(LDAfit@beta)
colnames(mybeta)<-LDAfit@terms
mybeta<-t(mybeta)
colnames(mybeta)<-seq(1:ncol(mybeta))
mybeta=exp(mybeta)

### First we print top 50 words
nwords=50
topwords <- mybeta[1:nwords,]
for (i in 1:LDAfit@k) {
  tempframe <- mybeta[order(-mybeta[,i]),]
  tempframe <- tempframe[1:nwords,]
  tempvec<-as.vector(rownames(tempframe))
  topwords[,i]<-tempvec
}
rownames(topwords)<-c(1:nwords)
write.xlsx(topwords, "TopWords.xlsx")

### Print top 30 documents
metadf<-data3
# notice that the "text" column is again named "text". If column name is different, name "text" needs to be changed.
meta_theta_df<-cbind(metadf[,"text"],LDAfit@gamma)
ntext=30
toptexts <- mybeta[1:ntext,]
for (i in 1:LDAfit@k) {
  print(i)
  tempframe <- meta_theta_df[order(-as.numeric(meta_theta_df[,i+1])),]
  tempframe <- tempframe[1:ntext,]
  tempvec<-as.vector(tempframe[,1])
  toptexts[,i]<-tempvec
}
rownames(toptexts)<-c(1:ntext)
write.xlsx(toptexts, "TopTexts.xlsx")


### Extrating unique words for topic (FREX words)
mybeta<-data.frame(LDAfit@beta)
colnames(mybeta)<-LDAfit@terms
mybeta<-t(mybeta)
colnames(mybeta)<-seq(1:ncol(mybeta))
mybeta=exp(mybeta)

# change myw to change the weight given to uniqueness
myw=0.3
word_beta_sums<-rowSums(mybeta)
my_beta_for_frex<-mybeta
for (m in 1:ncol(my_beta_for_frex)) {
  for (n in 1:nrow(my_beta_for_frex)) {
    my_beta_for_frex[n,m]<-1/(myw/(my_beta_for_frex[n,m]/word_beta_sums[n])+((1-myw)/my_beta_for_frex[n,m]))
  }
  print (m)
}
nwords=50
topfrex <- my_beta_for_frex[1:nwords,]
for (i in 1:LDAfit@k) {
  tempframe <- my_beta_for_frex[order(-my_beta_for_frex[,i]),]
  tempframe <- tempframe[1:nwords,]
  tempvec<-as.vector(rownames(tempframe))
  topfrex[,i]<-tempvec
}
rownames(topfrex)<-c(1:nwords)
write.xlsx(topfrex, "TopFREXWords.xlsx")


# Creating the network:
### names form data interpretation based on excel files
mynames<-c('Clinton Emails','Trump Bus Tape','Presidential Polls','Misc and Orlando','Artifact ','Artifact ','Economy and Education',
           'Campaign Data','Domestic Issues','GAO I','Artifact ','Artifact ','Horse Race','Ecigs','Presidential Race','Opioid Epidemic',
           'GOP Convention','Legislation Releases','Veterans Affairs','Artifact','Campaign Financing','PACs','Trump vs Kahn','Johnson Vs Feingold',
           'GAO II','Gun Control Legislation','Artifact  ','Unions Support','National Security','Foreign Affairs','Primaries','FEMA',
           'The Trump Nomination','SmartCity Grant','Opioid Legislation','Drug Enforcement','Artifact','Agency Regulation','Voting','Funding',
           'Trump vs GOP','Engineering Jobs','Obamacare','George Voinovich','Orlando Mass Shooting','Ohio Opioids','Battle States','Committee Reports',
           'New Hampshire Race','Obama','Police and crime','Labor Department','Mike Pence','Debates','Keep Guns from Terrorists',
           'Democrats Senate Hopes','Foreign Trade','Companies ','Climate Change','Paul Ryan and Trump','Medicare','Supreme Court',
           'Presidential campaign','Voters','Water regulation','Congress activity')

### using the network from LDA function:
mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           deleted_topics=c(5,6,11,12,20,27,37),
                           topic_names=mynames,
                           save_filename="trythis")

# We can also add the size of topics to the node attribute. In our example to improve model quality we removed duplicate entries
# however, we want to re-introduce these duplicates when calculating the topic salience. 
# If two candidates, for example, write the same message we want to include both their messages in topic salience, despite us removing it previously
# This might not be the case if duplicates are the result of error in data retrieval.
# if there is no need to re-populate duplicate documents the following code can be used:
LDAfit<-LDA.66
dfm_forsize<-data.frame(dfm_counts2)
dfm_forsize<-dfm_forsize[,-1]
sizevect<-rowSums(dfm_forsize)
meta_theta_df<-data.frame(size=sizevect,LDAfit@gamma)

topic.frequency <- colSums(meta_theta_df[,2:ncol(meta_theta_df)]*as.vector(meta_theta_df[,1]))
topic.proportion <- topic.frequency/sum(topic.frequency)

mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           deleted_topics=c(5,6,11,12,20,27,37),
                           topic_names=mynames,
                           save_filename="trythis",
                           topic_size = topic.proportion)

##### HOWEVER, if re-populating of duplicate topics is needed use the following code:
# we will use the theta data (topic*document matrix) from existing documents to assess the duplicated documents previously removed 
# first we prepare the meta data of existing documents and calculate their essential word count (words included in the topic model vocabulary after preprocessing and trimming)
LDAfit<-LDA.66
metadf<-data3
meta_theta_df<-cbind(metadf,LDAfit@gamma)
dfm_forsize<-data.frame(dfm_counts2)
dfm_forsize<-dfm_forsize[,-1]
sizevect<-rowSums(dfm_forsize)
meta_theta_df<-data.frame(size=sizevect,meta_theta_df)

# now we prepare the removed duplicates dataset
duplicate_df<-removed_df
colnames(duplicate_df)<-paste0(colnames(duplicate_df),".1")

# we cycle through all removed documents to add the missing theta values
dflist<-list()
for (i in (1:nrow(duplicate_df))) {
  the_match<-match(duplicate_df$text.1[i],meta_theta_df$text)
  newvect<-c(duplicate_df[i,],meta_theta_df[the_match,])
  dflist[[i]]<-newvect
}
maintable<-data.frame(do.call(bind_rows,dflist))

# we now delete the metadata from orginal matched document - leaving only meta data for the actual document with the theta values and size 
maintable<-data.frame(size=maintable$size,maintable[,-c((ncol(duplicate_df)+1):(ncol(duplicate_df)+ncol(metadf)+1))])
colnames(maintable)<-gsub("\\.1","",colnames(maintable))
meta_theta_df<-bind_rows(meta_theta_df,maintable)

# recalculating topic size
topic.frequency <- colSums(meta_theta_df[,(ncol(duplicate_df)+2):ncol(meta_theta_df)]*as.vector(meta_theta_df[,1]))
topic.proportion <- topic.frequency/sum(topic.frequency)

# using the function:
mynewnet<-network_from_LDA(LDAobject=LDAfit,
                           deleted_topics=c(5,6,11,12,20,27,37),
                           topic_names=mynames,
                           save_filename="trythis",
                           topic_size = topic.proportion)



