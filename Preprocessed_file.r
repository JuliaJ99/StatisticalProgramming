setwd("C:\\Users\\julia\\Desktop\\Julia_Git_repo\\StatisticalProgramming")
  
start_time= Sys.time()
a <- scan("pg10.txt",what="character",skip=104)
n <- length(a)
a <- a[-((n-2886):n)]
a <- a[-grep("[0123456789]:[0123456789]",a)]

# Separate out punctuation
split_punct <- function(x) {
  sub('(.)([,.;!:?])','\\1£\\2',x) |> 
    paste('£',sep="") |> 
    paste(collapse="") |>
    strsplit('£')
}

a <- split_punct(a)[[1]]
a[0:30]

# ------- #

# (a) # Determine unique words

unique_words = unique(tolower(a))
unique_words[0:10]

#(b) # The index of each word in a w.r.t unique_words
index_vector = match(tolower(a),unique_words)
index_vector

#(c) # Frequency of each word
tab = tabulate(index_vector)
tab

#(d) # Most frequent words
head(sort(tab,decreasing=TRUE),n=500)
###the word that occurs 163 times is the least of the 500 most common words.

#(e) # List of top 500 words
b = (unique_words[tab>=163])
b 

# -------- Task 7: Generating matrix T -------------------

# (a)
index_vector_freq=match(tolower(a),b)
index_vector_freq[0:20]

vec2 = index_vector_freq[2:length(index_vector_freq)] |> append(NA)
vec3 = index_vector_freq[3:length(index_vector_freq)] |> append(NA) |> append(NA)

col_matrix=cbind(index_vector_freq,vec2,vec3)

col_matrix_T= col_matrix[is.na(rowSums(col_matrix))== FALSE,1:3]

# b may contain words dropped after NA, just using it for now
T_matrix = array(0,dim=c(500, 500,500))
T_matrix

for (row in 1:length(col_matrix_T[,1])){
  current_row=col_matrix_T[row,]
  i=current_row[1]
  j=current_row[2]
  k=current_row[3]
  T_matrix[i,j,k]=T_matrix[i,j,k]+1
}

T_matrix[105,2,107]

# --------- Creating matrix A (for pairs of words)---------

col_matrix_A= cbind(index_vector_freq,vec2)
col_matrix_A

col_matrix_A= col_matrix_A[is.na(rowSums(col_matrix_A))== FALSE,]

A_matrix = array(0,dim=c(500,500)) #Should it be 500 or something else?
A_matrix

for (row in 1:length(col_matrix_A[,1])){
  current_row=col_matrix_A[row,]
  i=current_row[1]
  j=current_row[2]
  A_matrix[i,j]=A_matrix[i,j]+1
}

# ------------ Creating vector S
vector_S = tabulate(index_vector_freq)

# --------------Simulation
end_time=Sys.time()
totaltime=end_time-start_time

totaltime
