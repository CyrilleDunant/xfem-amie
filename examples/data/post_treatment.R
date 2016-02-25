######################
# Geometry functions #
######################

# computes the distance between two points
distance <- function( x1, y1, x2, y2 )
{
     return(sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )) ;
}

# computes the area of a triangle defined by the coordinates of its three nodes
triangle_area <- function( pts )
{
    a = distance(pts[1],pts[2],pts[3],pts[4]) ;
    b = distance(pts[3],pts[4],pts[5],pts[6]) ;
    c = distance(pts[5],pts[6],pts[1],pts[2]) ;
    s = (a+b+c)/2 
    return(sqrt(s*(s-a)*(s-b)*(s-c))) ;
}

# returns [ xmin, xmax, ymin, ymax ]
get_bounding_box <- function(data)
{
    return(c(min( data[,c(1,3,5)] ), max( data[,c(1,3,5)] ), min( data[,c(2,4,6)] ), max( data[,c(2,4,6)] ))) ;
}

# returns the dot-product between two points 
dot_product <- function( a, b )
{
    return( a[1]*b[1]+a[2]*b[2] ) ;
}

# returns the cross-product between two points 
cross_product <- function( a, b )
{
    return( a[1]*b[2]-a[2]*b[1] ) ;
}

# returns the intersection between two segments
segment_intersection <- function( seg_i, seg_j )
{
    inter=c() ;

    p = seg_i[1:2] ;
    r = seg_i[3:4]-p ;
    q = seg_j[1:2] ;
    s = seg_j[3:4]-q ;

    if(cross_product(r,s) == 0 & cross_product(q-p, r) == 0)
    {
        t0 = dot_product(q-p, r)/dot_product(r,r) ;
        t1 = t0 + dot_product(s,r)/dot_product(r,r) ;
        if(t0 >= 0 & t0 <= 1)
        {
            inter = c(inter, p[1]+t0*r[1], p[2]+t0*r[2]) ;
        }
        if(t1 >= 0 & t1 <= 1)
        {
            inter = c(inter, p[1]+t1*r[1], p[2]+t1*r[2]) ;
        }
        if(min(t0,t1) < 0 & max(t0,t1) >= 0)
        {
            inter = c(inter, p[1], p[2]) ;
        }
        if(min(t0,t1) <= 1 & max(t0,t1) > 1)
        {
            inter = c(inter, p[1]+r[1], p[2]+r[2]) ;
        }
    }
    if(cross_product(r,s) != 0)
    {
        t = cross_product(q-p,s)/cross_product(r,s) ;
        u = cross_product(q-p,r)/cross_product(r,s) ;
        if(t >= 0 & t <= 1 & u >= 0 & u <= 1)
        {
             inter = c(p[1]+t*r[1], p[2]+t*r[2]) ;
        }
        
    }
    return(inter) ;
}

# returns the intersection points bewteen a segment and a triangle
segment_triangle_intersection <- function( triangle, segment )
{
    inter_1 = segment_intersection( triangle[c(1,2,3,4)], segment ) ;
    inter_2 = segment_intersection( triangle[c(3,4,5,6)], segment ) ;
    inter_3 = segment_intersection( triangle[c(5,6,1,2)], segment ) ;

    if(length(inter_1[1]) == 4) { return(inter_1) ; }
    if(length(inter_2[1]) == 4) { return(inter_2) ; }
    if(length(inter_3[1]) == 4) { return(inter_3) ; }
    if(length(inter_1)+length(inter_2)+length(inter_3) == 6)
    {
        inter=c(inter_1[1], inter_1[2], inter_2[1], inter_2[2])
        if(inter_1[1]==inter_2[1] & inter_1[2]==inter_2[2])
        {
            inter[3] = inter_3[1] ;
            inter[4] = inter_3[2] ;
        }
        return(inter) ;
    }
    inter = inter_1 ;
    if(length(inter_2) == 2)
    {
        if(length(inter) == 0) { inter = inter_2 ; }
        if(length(inter) == 2 & (inter_2[1] != inter[1] | inter_2[2] != inter[2]) ) { inter = c(inter, inter_2[1], inter_2[2] ) ; }
    }
    if(length(inter_3) == 2)
    {
        if(length(inter) == 0) { inter = inter_3 ; }
        if(length(inter) == 2 & (inter_3[1] != inter[1] | inter_3[2] != inter[2]) ) { inter = c(inter, inter_3[1], inter_3[2] ) ; }
    }
    return(inter) ;
}

#######################
# Vector manipulation #
#######################

# returns a vector containing the middle points of a list of intervals
mid_interval <- function( intervals )
{
    mid = seq(1:(length(intervals)-1))*0 ;
    for(i in 1:length(mid))
    {
        mid[i] = (intervals[i]+intervals[i+1])/2
    }
    return(mid) 
}

# sorts the values of [y] by corresponding increasing values of [x]
sort_y <- function(x, y)
{
    sx = sort(x) ;
    sy = sx*0 ;
    for(i in 1:length(sx))
    {
        sy[i] = y[ x == sx[i] ] ;
    }
    return(sy) ;
}


# checks if a vector contains a certain value
contains <- function( vector, value )
{
    return( length(which(vector == value)) > 0 ) ;
}

# utility, compares two numbers using [operator]
compare <- function(a,b, operator="==")
{
    if(operator==">") { return (a > b) ; }
    if(operator==">=") { return (a >= b) ; }
    if(operator=="<") { return (a < b) ; }
    if(operator=="<=") { return (a <= b) ; }
    return(a == b) ;
}

################# 
# Index finders #
#################

# gets the index of the lines in [data] which correspond to triangles with [target] nodes located at a certain [y] coordinate
index_at_y <- function(data, y, operator="==", target=2)
{
    index = c() ;
    for( i in 1:dim(data)[1] )
    {
        count = 0 ;
        if(compare(data[i,2],y,operator)) { count = count+1 ; }
        if(compare(data[i,4],y,operator)) { count = count+1 ; }
        if(compare(data[i,6],y,operator)) { count = count+1 ; }
        if(count >= target)
        {
            index = c(index, i) ;
        }
    }
    return(index) ;
}

# gets the index of the lines in [data] which correspond to triangles with [target] nodes located at a certain [x] coordinate
index_at_x <- function(data, x, operator="==", target=2)
{
    index = c() ;
    for( i in 1:dim(data)[1] )
    {
        count = 0 ;
        if(compare(data[i,1],x,operator)) { count = count+1 ; }
        if(compare(data[i,3],x,operator)) { count = count+1 ; }
        if(compare(data[i,5],x,operator)) { count = count+1 ; }
        if(count >= target)
        {
            index = c(index, i) ;
        }
    }
    return(index) ;
}

# gets the index of the lines in [data] which correspond to triangles with the value of column [c] compared to [val]
index_at_c <- function(data, c, val, operator="==")
{
    index = c() ;
    for( i in 1:dim(data)[1] )
    {
        if(compare(data[i,c],val,operator))
        {
            index = c(index, i) ;
        }
    }
    return(index) ;
}

# gets the index of all lines in [data] which correspond to triangles that have at least [inter] intersection points with a [segment]
index_on_segment <- function( data, segment, inter=2)
{
    index = c() ;
    trg=seq(1,6) ;
    for( i in 1:dim(data)[1] )
    {
        if( length(triangle_segment_intersection( data[i,trg], segment )) >= inter*2 )
        {
            index = c(index, i) ;
        }
    }
    return(index) ;    
}

# returns all line from the data
index_all <- function(data)
{
    return( seq(1:dim(data)[1]) ) ;
}

# gets the union between two index lists (i OR j)
index_union <- function(index_i, index_j)
{
    return( unique(c(index_i, index_j))) ;
}

# gets the intersection between two index lists (i AND j)
index_intersection <- function(index_i, index_j)
{
    all = index_union(index_i, index_j) ;
    index = c() ;
    for(i in 1:length(all))
    {
        if( contains(index_i, all[i]) & contains(index_j, all[i]) )
        {
            index = c(index, all[i]) ;
        }
    }
    return(index) ;
}

# gets the complementary of [index] in [all]
index_complementary <- function( index, all )
{
    idx = c() ;
    for(i in 1:length(all))
    {
        if( length( which(index==all[i]) ) == 0 ) { idx = c(idx, all[i] ) ; }
    }
    return(idx) ;
}


##################
# Extract values #
##################

# returns the values of a column in data for all triangles which correspond to the selected index
# if multiple columns are called, their values will be averaged
values_at_index <- function(data, index, columns)
{
    values = seq(1:length(index))*0 ;
    for(i in 1:length(index))
    {
         for(j in 1:length(columns) )
         {
             values[i] = values[i] + (data[index[i], columns[j]]/length(columns)) ;
         }
    }
    return(values) ;
}

# returns the maximum distance between the center of a triangle pointed to by [index] and the selected [edge]
max_distance_from_edge <- function( data, index, edge="left")
{
    bbox = get_bounding_box(data) ;
    distance = -1 ;
    for(i in 1:length(index))
    {
        j = index[i] ;
        center = (data[j,1]+data[j,3]+data[j,5])/3 ;
        if( edge=="top" | edge=="bottom") { center = (data[j,2]+data[j,4]+data[j,6])/3 ; }
        current = -1 ;
        if( edge == "left" ) { current = center-bbox[1] ; }
        if( edge == "right" ) { current = bbox[2]-center ; }
        if( edge == "bottom" ) { current = center-bbox[3] ; }
        if( edge == "top" ) { current = bbox[4]-center ; }
        if(current > distance) { distance = current ; }
    }
    return(distance) ;
}

##########################
# Extract average values #
##########################

# returns the surface average of [column] over all triangles pointed to by [index]
surface_average <- function(data, index, column)
{
    values = 0 ;
    weights= 0 ;

    for(i in 1:length(index))
    {
        j = index[i] ;
        v = data[j,column] ;
        a = triangle_area( c(data[j,1],data[j,2],data[j,3],data[j,4],data[j,5],data[j,6]) ) ;
        
        values = values+(v*a) ;
        weights = weights+a ;
    }

    if(weights > 0)
    {
        values = values / weights ;
    }
    return(values) ;
    
}

# returns the surface average of [column] over a series of [slices] taken along a certain [direction]
profile_surface_average <- function( data, slices, column, direction="x")
{
    values = array(0,length(slices)-1) ;
    weights = array(0,length(slices)-1) ;

    for(i in 1:dim(data)[1])
    {
        v = data[i,column] ;
        a = triangle_area( c(data[i,1],data[i,2],data[i,3],data[i,4],data[i,5],data[i,6]) ) ;
        center = (data[i,1]+data[i,3]+data[i,5])/3 ;
        if(direction == "y")
        {
            center = (data[i,2]+data[i,4]+data[i,6])/3 ;
        }
        for( j in 1:length(values) )
        {
            if( center >= slices[j] & center <= slices[j+1] )
            {
                values[j] = values[j] + v*a ;
                weights[j] = weights[j] + a ;
                break ;
            }
        }
    }

    for(i in 1:length(values))
    {
        if(weights[i] > 0)
        {
            values[i] = values[i]/weights[i] ;
        }
    }
    return(values) ;
}

# returns the surface average of [column] over a series of [slices] taken along a certain [direction], accounting only for elements in [index]
profile_surface_average_in_index <- function( data, index, slices, column, direction="x")
{
    values = array(0,length(slices)-1) ;
    weights = array(0,length(slices)-1) ;

    for(ii in 1:length(index) )
    {
        i = index[ii] ;
        v = data[i,column] ;
        a = triangle_area( c(data[i,1],data[i,2],data[i,3],data[i,4],data[i,5],data[i,6]) ) ;
        center = (data[i,1]+data[i,3]+data[i,5])/3 ;
        if(direction == "y")
        {
            center = (data[i,2]+data[i,4]+data[i,6])/3 ;
        }
        for( j in 1:length(values) )
        {
            if( center >= slices[j] & center <= slices[j+1] )
            {
                values[j] = values[j] + v*a ;
                weights[j] = weights[j] + a ;
                break ;
            }
        }
    }

    for(i in 1:length(values))
    {
        if(weights[i] > 0)
        {
            values[i] = values[i]/weights[i] ;
        }
    }
    return(values) ;
}

#############################
# Post-treatment operations #
#############################

# returns the total area of the triangles pointed to by [index]
get_index_area <- function(data, index=c())
{
    surface = index*0 ;
    trg=seq(1,6) ;
    if(length(surface) == 0) { surface = array(0, dim(data)[1]) ; }
    for(i in 1:length(surface))
    {
        j = i ;
        if( length(index) > 0) { j = index[i] ; }
        surface[i] = triangle_area( data[j,trg] ) 
    }
    return(surface) ;
}

# extract principal components from 3 vectors
get_principal_values <- function(v11, v22, v12)
{
    principal = array(0,c(length(v11), 2)) ;

    for(i in 1:length(v11))
    {
        trace = v11[i] + v22[i] ;
        det = v11[i]*v22[i] - 0.25*v12[i]*v12[i] ;
        delta = sqrt(trace*trace-4*det) ;
        angle = 0.5*atan2( v11[i]-v22[i], 0.5*v12[i] ) ;
        if(cos(angle) < 0)
        {
            principal[i,1] = (trace + delta)*0.5 ;
            principal[i,2] = (trace - delta)*0.5 ;
        }
        else
        {
            principal[i,1] = (trace - delta)*0.5 ;
            principal[i,2] = (trace + delta)*0.5 ;
        }
    }
    return(principal) ;
} 

