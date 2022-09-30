setwd("C:\\Users\\julia\\Desktop\\Julia_Git_repo\\StatisticalProgramming")
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


