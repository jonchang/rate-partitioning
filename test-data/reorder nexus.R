#transpose data
data.t <- t(data)

#sort by column 1 (rel rate)
data.t[rev(order(data.t[, 1])), ]

# makes data frame look better by getting ride of quote marks?
df <- as.data.frame(data.t[rev(order(data.t[, 1])), ])

#transposes data back into column as nucleotide sites
t(df)

#cluster data into k bins where k=2
kmeans(as.numeric(as.character(df$V1)), 2)
