import numpy as np

def lookup(lookupval,lookuparray,resultarray,threshold):
    #defining a lookup function that assumes lookuparray is sorted in ascending 
    #order
    
    #**threshold is the maximum difference between the lookupval and the value
    #found in the lookup aray in order for the function to return a result. For
    #example, in abundance matching, we would typically set threshold to the
    #average step in numden.
    
    lookuparval,i=float('-inf'),0
    #while the values in the lookup array are less than the lookup value:  
    while lookuparval<lookupval: 
        if i>lookuparray.size-1: 
            #If the loop runs through every element of the lookuparray and 
            #lookuparval is still < lookupval, end the loop, otherwise an error 
            #will occur:
            i-=1
            #Work backwards to find i corresponding to the last time 
            #lookuparval increased:
            while lookuparray[i]==lookuparray[i-1]:
                i-=1
            #Once the nested loop has gotten i back to the point where 
            #lookuparray[i]!=lookuparray[i], break the primary loop:
            break 
        lookuparval=lookuparray[i]
        
        i+=1
        if i>1e4:
            sys.exit("Timeout")
    lowerlookuparval=lookuparray[max(0,i-2)]
    upperlookuparval=lookuparray[max(0,i-1)]
    
    #Check whether the lookuparval below the lookupval or the one above is 
    #closer to the lookupval and use the one that's closer
    upperdiff=np.abs(upperlookuparval-lookupval)
    lowerdiff=np.abs(lookupval-lowerlookuparval)
    if upperdiff>lowerdiff:
        #If lookuparval below the lookupval is closer, make sure lowerdiff is
        #within the threshold.
        if np.abs(lowerdiff)<threshold: 
            result=resultarray[max(0,i-2)]
            print('halo number density: {0:0.2f}'.format(lowerlookupval))
        else:
            #return #EXCLUDING "NA" FOR NOW
            result=None
    #If lookuparval above the lookupval is closer, make sure upperdiff is 
    #within the threshold.
    elif np.abs(upperdiff)<threshold:
        result=resultarray[max(0,i-1)]
        print('halo number density: {0:0.2f}'.format(upperlookupval))
    else:
        #return #EXCLUDING "NA" FOR NOW
        result=None
    
    return result
