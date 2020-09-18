# fast matrix products using C code
mv <- function(A, B, transpose=FALSE){
    if(!transpose){
        as.matrix(.Call("mv_c",  A, B))
    } else if(transpose){
        as.matrix(.Call("tmv_c", A, B))
    }
}
