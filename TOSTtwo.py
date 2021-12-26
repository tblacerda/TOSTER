
def TOSTtwo(m1, m2,
            sd1, sd2,
            n1, n2,
            low_eqbound_d,
            high_eqbound_d,
            alpha = 0.05,
            var_equal = False,
            plot = False,
            verbose = True):

    if (n1 < 2) or (n2 < 2):
        return "The sample size should be larger than 1."

    if (1<=alpha or alpha < 0):
        return "The alpha level should be a positive value between 0 and 1."
    
    if (sd1 <= 0 or sd2 <=0):
        return "The standard deviation should be a positive value."
    
    ## Fim dos checks
      # Calculate TOST, t-test, 90% CIs and 95% CIs
      
    if var_equal == True:
        sdpooled = sqrt((((n1 - 1)*(sd1**2))+(n2 - 1)*(sd2**2))/((n1+n2)-2))
        low_eqbound = low_eqbound_d*sdpooled
        high_eqbound = high_eqbound_d*sdpooled
        degree_f = n1+n2-2
    
        dist = student_t(df=degree_f,loc=0,scale=1 )


        t1 = ((m1-m2)-low_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))  #students t-test lower bound
        lower_tail_false = 1- dist.cdf(t1)  
        p1 = lower_tail_false 
        t2 = ((m1-m2)-high_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2)) #students t-test upper bound
        lower_tail_true = dist.cdf(t2)
        p2 = lower_tail_true
        
        t = (m1-m2)/(sdpooled*sqrt(1/n1 + 1/n2))
        
        lower_tail_true2 = dist.cdf(-abs(t))
        pttest = 2*lower_tail_true2
        
        LL90 = (m1-m2)-student_t.ppf(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
        UL90 = (m1-m2)+student_t.ppf(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
        LL95 = (m1-m2)-student_t.ppf(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
        UL95 = (m1-m2)+student_t.ppf(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    else:
        sdpooled = sqrt((sd1**2 + sd2**2)/2) #calculate sd root mean squared for Welch's t-test
        low_eqbound = low_eqbound_d*sdpooled
        high_eqbound = high_eqbound_d*sdpooled
        degree_f = (sd1**2/n1+sd2**2/n2)**2/(((sd1**2/n1)**2/(n1-1))+((sd2**2/n2)**2/(n2-1))) #degrees of freedom for Welch's t-test        
        dist = student_t(df=degree_f,loc=0,scale=1 )
        t1 = ((m1-m2)-low_eqbound)/sqrt(sd1**2/n1 + sd2**2/n2) #welch's t-test upper bound
        lower_tail_false = 1- dist.cdf(t1)  
        p1 = lower_tail_false 
        t2 = ((m1-m2)-high_eqbound)/sqrt(sd1**2/n1 + sd2**2/n2) #welch's t-test lower bound
        lower_tail_true = dist.cdf(t2)
        p2 = lower_tail_true
        t = (m1-m2)/sqrt(sd1**2/n1 + sd2**2/n2) #welch's t-test NHST    
        lower_tail_true2 = dist.cdf(-abs(t))
        pttest = 2*lower_tail_true2
    
        LL90 = (m1-m2)-student_t.ppf(1-alpha, degree_f)*sqrt(sd1**2/n1 + sd2**2/n2) #Lower limit for CI Welch's t-test
        UL90 = (m1-m2)+student_t.ppf(1-alpha, degree_f)*sqrt(sd1**2/n1 + sd2**2/n2) #Upper limit for CI Welch's t-test
        LL95 = (m1-m2)-student_t.ppf(1-(alpha/2), degree_f)*sqrt(sd1**2/n1 + sd2**2/n2) #Lower limit for CI Welch's t-test
        UL95 = (m1-m2)+student_t.ppf(1-(alpha/2), degree_f)*sqrt(sd1**2/n1 + sd2**2/n2) #Upper limit for CI Welch's t-test
  
    ptost = max(p1,p2) #Get highest p-value for summary TOST result
    ttost = t2
    if (abs(t1) < abs(t2)):
        ttost = t1
  
    dif = (m1-m2)
    testoutcome = "non-significant"
    
    if pttest < alpha:
        testoutcome = "significant"
    
    TOSToutcome = "non-significant"
    if ptost<alpha:
        TOSToutcome = "significant"
    
    print("p1= ", p1)
    print("p2= ", p2)
    print("t1= ", t1)
    print("t2= ", t2)
    print("ptost= ", ptost)
    print("ttost= ", ttost)
    print('testOutcome = ', testoutcome)
    print('TOSToutcome = ', TOSToutcome)
    print('Lower Conf. Interval = ', LL90)
    print('Upper Conf. Interval = ', UL90)
    
    return np.nan

TOSTtwo(m1 = 4.55,
        m2 = 4.87,
        sd1 = 1.05,
        sd2 = 1.11,
        n1 = 150,
        n2 = 15,
        low_eqbound_d= -0.5,
        high_eqbound_d= 0.5)

