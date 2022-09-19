rm(list=ls())
## Do Ctrl Enter to run each line 
#### Arithmetic ####
17 + 3 ## Addition
17 - 3 ##Subtraction
17 * 3 ##Multiplication
17 / 3 ##Division
17 ^ 3 ##Exponentiation
17 %% 3 ##Modulo (remainder from division)
17 %/% 3 ##Integer Division 

##### Equalities ####
17 == 3 ##Equal to
17 != 3 ##Not equal to
17 > 3 ##Greater than
3 < 17 ##Less than
17 >= 3 ##Greater than or equal 
17 <= 3 ##Less than or equal
17 > 3 | 17 < 3 ##| is or
17 > 3 & 17 < 3 ##& is and 

##### Variables ####
##clear all variables by rm(list=ls())
x <- -1.2345
y <- "Hi"
z <- 6.34 * 10^5
x
y
z
rm(list=ls())

##### Data Types ####
##INteger: whole numbers
##numeric: real numbers
##complex: complex numbers
##logical: true or false
##character: string
##specifiy data type using as.type(variable)
##query data type use is.type(variable)
x <- as.integer(50) 
class(x)
1+x
z <- as.numeric(x^2)
y <- as.character("Bye")
q <- as.complex(1+3i)
r <- 3==4
is.numeric(x)
is.numeric(y)
is.numeric(z)
is.numeric(q)
is.numeric(r)
is.logical(r)

#### Math Functions ####
rm(list=ls())
x <- -1.234
abs_x <- abs(x) #absolute value
sqrt_abs_x <- sqrt(abs_x) #square root
ceiling(x); ceiling(abs_x) #round up
floor(x); floor(abs_x) #round down
trunc(x); trunc(abs_x)#remove decimals
(y <- round(x, 2))
(cos_y <- cos(y))
(z <- log(abs(y)))
(exp_z <- exp(z))

#### Vectors is a collection of elements####
rm(list=ls())
1:6
x <- c(1:6)
y <- c("a", "B", "c", "D")
y[3] #asks for third position in vector y
z <- seq(1, 40, 4) #first one is where to start the sequence and second one is where to end the sequence; third one is adding 4
z
z[3:6]
z[c(2,8,10)]

#### Statistical functions for vectors ####
length(y) # number of elements
min(x); max(x)
sum(x); prod(x)
median(z)
mean(z)
var(z)
sd(z) #standard deviation
IQR(z) #inner quartile ranges 
summary(z)



#Base R Mind Expander Exercise
my_vector <- seq(2,100,2) ##makes a sequence of 2 to 100 for even numbers 

##shows contents of my vectors
my_vector
(my_vector <- seq(2,100,2)) #shows contents of my_vector
print(my_vector)
view(my_vector)
summary(my_vector)

#showing all elements in my_vector divisible by 12
my_vector%%12 == 0
my_vector_div_12 <- my_vector[my_vector%%12 == 0] ##shows all elements divisible by 12 by moduloing the ones (having a remainder) in my_vector that have 0 remainder
my_vector_div_12
sum(my_vector_div_12) ##shows sum of all the elements in my_vector_div_12

prod(my_vector[c(5,10,15)]) ##extracts elements 5,10,15 of my_vector and multiplies all of them

