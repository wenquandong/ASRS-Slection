#!/usr/bin/env python
# coding: utf-8

# In[1]:


def threshholds(NA):
    threshhold=[max(Q.values())]
    Cst={k:min(C_s*n+C_v+C_t*n*k*min(p[i]/Q[i] for i in Q)*(EC31[k]),
               C_s+C_v+C_t*n*k*min(p[i]/Q[i] for i in Q)*(1/n*EC31[k]+(1-1/n)*EC32[k]))for k in range(2,max(Q.values())+1)}
    Cst[1]=C_v+C_t*n*1*min(p[i]/Q[i] for i in Q)*(EC31[1])
    for k in range(3,1+max(Q.values())):
        a1=2
        a2=k-2
        b1=math.floor(k/2)
        b2=math.ceil(k/2)
        if Cst[k]>=Cst[b1]+Cst[b2]:
            threshhold.append(k)
    remove=[]
    for k in threshhold:
        if sum(math.ceil(Q[i]/k) for i in Q)>n*NA:
            remove.append(k)
    #print (remove)
    for k in remove:
        threshhold.remove(k)
    threshhold.append(max(Q.values()))
    #K=int(min(threshhold))
    #print (K)
    return (threshhold)
def rootnode(NA):
    #startTime=datetime.now()
    m = Model("MODEL")
    y = {}
    n_ij = {}
    n_z = {}
    K=max(Q.values())
    a={}
    b={}
    C_st={}
    s={}
    for j in range(1,NA+1):
        a[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
        b[j]=m.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
        C_st[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
        s[j]=m.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
        for k in range(0,K+1):
            n_z[j,k]=m.addVar(vtype=GRB.BINARY, name = "n_z%s" % str([j,k]))
    for j in range(1,NA+1): 
        for i in range(1,I+1):
            y[i, j] = m.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
            n_ij[i,j]=m.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
    m.setObjective(quicksum(C_st[j]+C_v*quicksum(n_z[j,k] for k in range(1,K+1))
                                +C_rw*w*h*l*n*quicksum(k*n_z[j,k] for k in range(0,K+1)) for j in range (1,NA+1))
                       , GRB.MINIMIZE)

    for i in range(1,I+1):
        m.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
    for j in range(1,NA+1):
        m.addConstr(quicksum(n_z[j,k] for k in range(0,K+1))==1,"constriants 2")
    for j in range(1,NA+1):
        m.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
        for i in range(1,I+1):
            for k in range(0,K+1):
                m.addConstr(n_ij[i,j]*k>=y[i,j]-Q[i]*(1-n_z[j,k]), "constriants 11")
    for j in range(1,NA+1):
        for k in range(2,K+1):
            
            m.addConstr(a[j]>=C_s+(1/n*EC31[k]+(1-1/n)*EC32[k])*C_t*
                        quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))-M3*(1-n_z[j,k]) ,"new constriants 3")
            m.addConstr(b[j]>=C_s*n+EC31[k]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                        -M3*(1-n_z[j,k]),"new constriants 4")
            m.addConstr(a[j]-b[j]<=M3*s[j],"new constriants 5")
            m.addConstr(b[j]-a[j]<=M3*(1-s[j]),"new constriants 6")
        m.addConstr(C_st[j]>=b[j]-M3*(1-s[j]+n_z[j,1]+n_z[j,0]),"constriants 7")
        m.addConstr(C_st[j]>=a[j]-M3*(s[j]+n_z[j,1]+n_z[j,0]),"constriants 8")
        m.addConstr(C_st[j]>=quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))*4/3*Txy*C_t,"constriants 9")
    for j in range(1,NA):
        m.addConstr(quicksum(n_z[j,k]*k for k in range(0,K+1))<=
                        quicksum(n_z[j+1,k]*k for k in range(0,K+1)),"constriants 12")

    #m.params.LogToConsole=0
    if NA>1:
        m.params.SolutionLimit=10
        m.params.TimeLimit=60
    #model.getEnv().set(GRB_IntParam_OutputFlag, 0)
   # model.write("zainali.lp")
    m.update()
    m.optimize()
    #print ((datetime.now() - startTime).total_seconds())
    solution1={(i,j):y[i,j].x for i in range(1,I+1) for j in range(1,NA+1)}
    solution2={(i,j):n_ij[i,j].x for i in range(1,I+1) for j in range(1,NA+1)}
    solution3={j:sum(k*n_z[j,k].x for k in range(0,K+1)) for j in range(1,NA+1) }
    solution4={j:0 for j in range(1,NA+1)}
    for j in range(1,NA+1):
        if s[j].x==0:
            for k in range(2,1+K):
                if n_z[j,k].x==1:
                    solution4[j]=1
        else:
            for k in range(2,1+K):
                if n_z[j,k].x==1:
                    solution4[j]=n
    #print (solution3)
    return (m.ObjVal,solution1,solution2,solution3,solution4)
    
def LB(klist): ###
    #startTime=datetime.now()
    point=[i for i in klist]
    Qv={i:p[i]/Q[i] for i in Q}
    Q1={}
    while len(Qv)>0:
        match=max(Qv, key=Qv.get)
        Q1[match]=Q[match]
        del Qv[match]
    Y={(i,j):0 for i in Q1 for j in range(1,len(point)+1)}
    for j in range(1,len(point)+1):
        for i in Q1:
            if sum(Y[(SKU,j)] for SKU in Q1)<n*point[j-1]:
                Y[(i,j)]=min(Q1[i],n*point[j-1]-sum(Y[(SKU,j)] for SKU in Q1 ))
                Q1[i]-=Y[(i,j)]
    cost=0
    for rack in range(1,len(point)+1):
        if point[rack-1]==1:
            cost+=C_v+C_rw*w*h*l*n+C_t*EC31[int(point[rack-1])]*sum(Y[i,rack]/Q[i]*p[i] for i in Q)
        else:
            cost+=C_v+C_rw*w*h*l*n*point[rack-1]+min(C_s*n+C_t*sum(Y[i,rack]/Q[i]*p[i] for i in Q)*EC31[int(point[rack-1])],
                                                     C_s+C_t*sum(Y[i,rack]/Q[i]*p[i] for i in Q)*
                                                    (1/n*EC31[int(point[rack-1])]+(1-1/n)*EC32[int(point[rack-1])]))
    if sum(Q1.values())>0:
        if point[-1]>=2:
        
            cost+=C_v+C_rw*w*h*l*(sum(Q1.values()))+min(C_s*n+C_t*sum(Q1[i]/Q[i]*p[i] for i in Q)*EC31[int(point[-1])],
                                                             C_s+C_t*sum(Q1[i]/Q[i]*p[i] for i in Q)*
                                                            (1/n*EC31[int(point[-1])]+(1-1/n)*EC32[int(point[-1])]))
        else:
            cost+=C_v+C_rw*w*h*l*(sum(Q1.values()))+C_t*sum(Q1[i]/Q[i]*p[i] for i in Q)*EC31[int(point[-1])]
                                                             
    #print ('LB',klist,math.floor(cost))
    return math.floor(cost)



def UB(klist1):
    
    #startTime=datetime.now()
    #K=max(Q.values())
    point1={i:klist1[i-1] for i in range(1,len(klist1)+1)}
    #print (klist1,point1,n)
    #rint (klist)
    #ub=int(max(klist.values()))
    NA=len(point1)
    if n*sum(point1.values())<sum(Q.values()):
        #print ('UB',klist1,(datetime.now() - startTime).total_seconds())
        return (math.inf)
    else:
       
        kprime=math.floor(Txy*V_s/2/l)
        M1=sum(Q[i] for i in range(1,I+1))
        M2=max(w*h*sum(Q.values())/V_cx/V_cy,2*K*l/V_s)
        M3=C_s*n+(1/n*EC31[K]+(1-1/n)*EC32[K])*C_t
        #startTime=datetime.now()
        M = Model("MODEL")
        # Creating Variables
        y = {}
        n_ij = {}
        n_z = {}
        m={}
        a={}
        b={}
        C_st={}
        s={}
        for j in range(1,NA+1):
            a[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
            b[j]=M.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
            C_st[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
            s[j]=M.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
        for j in range(1,NA+1): 
            for i in range(1,I+1):
                y[i, j] = M.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
                n_ij[i,j]=M.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
        M.setObjective(quicksum(C_st[j] for j in point1)+C_v*len(point1)
                           +C_rw*w*h*l*n*sum(point1.values()), GRB.MINIMIZE)

        for i in range(1,I+1):
            M.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
        for j in range(1,NA+1):
            M.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
            for i in range(1,I+1):
                    M.addConstr(n_ij[i,j]*point1[j]>=y[i,j], "constriants 11")

        #M1
        for j in range(1,NA+1):
            if point1[j]>=2:
                
                M.addConstr(a[j]>=C_s+(1/n*EC31[point1[j]]+(1-1/n)*EC32[point1[j]])*C_t*
                            quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1)),"new constriants 3")
                M.addConstr(b[j]>=C_s*n+EC31[point1[j]]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                           ,"new constriants 4")
                M.addConstr(a[j]-b[j]<=M3*s[j],"new constriants 5")
                M.addConstr(b[j]-a[j]<=M3*(1-s[j]),"new constriants 6")
                M.addConstr(C_st[j]>=b[j]-M3*(1-s[j]),"constriants 7")
                M.addConstr(C_st[j]>=a[j]-M3*s[j],"constriants 8")
            else:
                M.addConstr(C_st[j]>=quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))*4/3*Txy*C_t,"constriants 9")
        M.params.LogToConsole=0
        #M.params.SolutionLimit=1
        M.params.TimeLimit=60
        #model.params.IntFeasTol=1e-9
        M.reset()
        M.update()
        M.optimize()

        if M.status==GRB.INFEASIBLE:
            #print ('UB',klist1,(datetime.now() - startTime).total_seconds())
            return (math.inf)
        else:
            print ('UB',klist1,M.ObjVal)
            return (M.ObjVal)
            

# consider all four senario
def DFS_4(graph, vb):
    #print (vb)
    global zub
    global zlb
    global lo
    stack = []
    order=[]
    seen = set()
    seen.add(vb)
    #lo=min(LO[i] for i in range(0,treenode[vb][0]+1))
    if sum(treenode[vb][1])<=lo:
        ZLB[vb]=LB(treenode[vb][1])
        if ZLB[vb]<=zub:
            ZUB[vb]=UB(treenode[vb][1])
            if ZUB[vb]<=zub:
                lo=sum(treenode[vb][1])
                #stack.append(vb)
                zub=ZUB[vb]
                
                BEST=treenode[vb][1]
                print (BEST,zub)
            elif ZUB[vb]==math.inf and len(treenode[vb][1])<lo:
                stack.append(vb)
            if ZLB[vb]>zlb:
                zlb=ZLB[vb]
        #stack.append(s)
            while(len(stack) > 0):
                vertex = stack.pop()
                order.append(vertex)
                #print (vertex)
                nodes = sorted(graph[vertex],reverse=True)
                #print (nodes)
                for ws in nodes:
                    #3lo=min(LO[i] for i in range(0,treenode[ws][0]+1))
                    if ws not in seen and lo>=sum(treenode[ws][1]):
                    #if ws not in seen:
                        #print (w)
                        seen.add(ws)
                        ZLB[ws]=LB(treenode[ws][1])
                        if ZLB[ws]<=zub:
                            #stack.append(w)
                            ZUB[ws]=UB(treenode[ws][1])
                            if ZUB[ws]<=zub:
                                #lo=sum(treenode[ws][1])
                                #lo=sum(treenode[ws][1])
                                lo=sum(treenode[ws][1])
                                #stack.append(ws)
                                zub=ZUB[ws]
                                BEST=treenode[ws][1]
                                #print (BEST,zub)
                            elif ZUB[ws]==math.inf and len(treenode[ws][1])<NAmax:
                                stack.append(ws)  
                            if ZLB[ws]>zlb:
                                zlb=ZLB[ws]
    #print ('point discovered',len(order))
    return order
# only consider the first three senario
def DFS_3(graph, vb):
    #print (vb)
    global zub
    global zlb
    global LO
    stack = []
    order=[]
    seen = set()
    seen.add(vb)
    lo=min(LO[i] for i in range(0,treenode[vb][0]+1))
    if sum(treenode[vb][1])<=lo:
        ZLB[vb]=LB(treenode[vb][1])
        if ZLB[vb]<=zub:
            ZUB[vb]=UB(treenode[vb][1])
            if ZUB[vb]<=zub:
                LO[treenode[vb][0]]=sum(treenode[vb][1])
                #stack.append(vb)
                zub=ZUB[vb]
                
                BEST=treenode[vb][1]
                print (BEST,zub)
            elif ZUB[vb]==math.inf and len(treenode[vb][1])<lo:
                stack.append(vb)
            if ZLB[vb]>zlb:
                zlb=ZLB[vb]
        #stack.append(s)
            while(len(stack) > 0):
                vertex = stack.pop()
                order.append(vertex)
                #print (vertex)
                nodes = sorted(graph[vertex],reverse=True)
                #print (nodes)
                for ws in nodes:
                    lo=min(LO[i] for i in range(0,treenode[ws][0]+1))
                    if ws not in seen and lo>=sum(treenode[ws][1]):
                    #if ws not in seen:
                        #print (w)
                        seen.add(ws)
                        ZLB[ws]=LB(treenode[ws][1])
                        if ZLB[ws]<=zub:
                            #stack.append(w)
                            ZUB[ws]=UB(treenode[ws][1])
                            if ZUB[ws]!=math.inf:
                                LO[treenode[ws][0]]=min(sum(treenode[ws][1]),lo)
                            if ZUB[ws]<=zub:
                                
                                #lo=sum(treenode[ws][1])
                                #lo=sum(treenode[ws][1])
                                #stack.append(ws)
                                zub=ZUB[ws]
                                BEST=treenode[ws][1]
                                print (BEST,zub)
                            elif ZUB[ws]==math.inf and len(treenode[ws][1])<NAmax:
                                stack.append(ws)  
                            if ZLB[ws]>zlb:
                                zlb=ZLB[ws]
    #print ('point discovered',len(order))
    return order
##without considering total depth
def DFS(graph, vb):
    #print (vb)
    global zub
    global zlb
    global LO
    stack = []
    order=[]
    seen = set()
    seen.add(vb)
    #lo=min(LO[i] for i in range(0,treenode[vb][0]+1))
    #if sum(treenode[vb][1])<=lo:
    ZLB[vb]=LB(treenode[vb][1])
    if ZLB[vb]<=zub:
        ZUB[vb]=UB(treenode[vb][1])
        if ZUB[vb]<=zub:
            #LO[treenode[vb][0]]=sum(treenode[vb][1])
            stack.append(vb)
            zub=ZUB[vb]

            BEST=treenode[vb][1]
            #print (BEST,zub)
        #elif ZUB[vb]==math.inf and len(treenode[vb][1])<NAmax:
        elif ZUB[vb]==math.inf and len(treenode[vb][1])<len(dep)-1:
            stack.append(vb)
        if ZLB[vb]>zlb:
            zlb=ZLB[vb]
        #stack.append(s)
        while(len(stack) > 0):
            vertex = stack.pop()
            order.append(vertex)
            #print (vertex)
            nodes = sorted(graph[vertex],reverse=True)
            #print (nodes)
            for ws in nodes:
                #lo=min(LO[i] for i in range(0,treenode[ws][0]+1))
                if ws not in seen:
                #if ws not in seen:
                    #print (w)
                    seen.add(ws)
                    ZLB[ws]=LB(treenode[ws][1])
                    if ZLB[ws]<=zub:
                        #stack.append(w)
                        ZUB[ws]=UB(treenode[ws][1])
                        if ZUB[ws]<=zub:
                            #LO[treenode[ws][0]]=sum(treenode[ws][1])
                            #lo=sum(treenode[ws][1])
                            #lo=sum(treenode[ws][1])
                            stack.append(ws)
                            zub=ZUB[ws]
                            BEST=treenode[ws][1]
                            print (BEST,zub)
                        elif ZUB[ws]==math.inf and len(treenode[ws][1])<NAmax:
                            stack.append(ws)  
                        if ZLB[ws]>zlb:
                            zlb=ZLB[ws]
    #print ('point discovered',len(order))
    return order

    
    
    
    


# In[2]:


def Objective_value(klist1):
    point1={i:klist1[i-1] for i in range(1,len(klist1)+1)}
    #startTime=datetime.now()
    NA=len(point1)
    M = Model("MODEL")
    # Creating Variables
    y = {}
    n_ij = {}
    n_z = {}
    m={}
    a={}
    b={}
    C_st={}
    s={}
    for j in range(1,NA+1):
        a[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
        b[j]=M.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
        C_st[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
        s[j]=M.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
    for j in range(1,NA+1): 
        for i in range(1,I+1):
            y[i, j] = M.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
            n_ij[i,j]=M.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
    M.setObjective(quicksum(C_st[j] for j in point1)+C_v*len(point1)+C_rw*w*h*l*n*sum(point1.values()), GRB.MINIMIZE)

    for i in range(1,I+1):
        M.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
    for j in range(1,NA+1):
        M.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
        for i in range(1,I+1):
                M.addConstr(n_ij[i,j]*point1[j]>=y[i,j], "constriants 11")

    #M1
    for j in range(1,NA+1):
        if point1[j]>=2:
            M.addConstr(a[j]>=C_s+(1/n*EC31[point1[j]]+(1-1/n)*EC32[point1[j]])*C_t*
                        quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1)),"new constriants 3")
            M.addConstr(b[j]>=C_s*n+EC31[point1[j]]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                       ,"new constriants 4")
            M.addConstr(a[j]-b[j]<=M3*s[j],"new constriants 5")
            M.addConstr(b[j]-a[j]<=M3*(1-s[j]),"new constriants 6")
            M.addConstr(C_st[j]>=b[j]-M3*(1-s[j]),"constriants 7")
            M.addConstr(C_st[j]>=a[j]-M3*s[j],"constriants 8")
        else:
            M.addConstr(C_st[j]>=quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))*4/3*Txy*C_t,"constriants 9")
    M.params.LogToConsole=0
    M.params.TimeLimit=60
    M.reset()
    M.update()
    M.optimize()
    print (M.ObjVal, point1)
    return (M.ObjVal, point1)
def find_depth(NAmin,NAmax):  
    status=['none']
    solution=['none']
    solutiondep=[math.inf]
    obj={NA:math.inf for NA in range(NAmin-1,NAmax+1)}
    for NA in range(NAmin,NAmax+1):
        Txy=math.sqrt(w*h*n/V_cx/V_cy)
        kprime=math.floor(Txy*V_s/2/l)
        K=max(Q.values())
        #print (kprime)
        #threshhold=[max(Q.values())/max(1,math.floor(n/I))]
        Tz={k:2*l*k/V_s for k in range(0,1+max(Q.values()))}
        EC2=4/3*Txy
        EC31={} #EC3Dj part 1 without m/n
        EC32={} #EC3Dj part 2
        EC31[0]=0
        EC32[0]=0
        EC31[1]=4/3*Txy
        EC32[1]=4/3*Txy
        for k in range(2,max(Q.values())+1):
            EC32[k]=9/5*Txy+k*l/V_s
            if Txy>=2*k*l/V_s:
                EC31[k]=4/3*Txy+(2*k*l/V_s/Txy)**3/12*Txy
            else:
                EC31[k]=2/3*Txy+2*k*l/V_s*(1/2+(Txy/(2*k*l/V_s))**2/4)
        M1=sum(Q[i] for i in range(1,I+1))
        M2=max(w*h*sum(Q.values())/V_cx/V_cy,2*K*l/V_s)
        M3=C_s*n+(1/n*EC31[K]+(1-1/n)*EC32[K])*C_t
        #####################################
        K=min(threshholds(NA)) 
        model = Model("MODEL")
        #print ('new model',NA)
        # Creating Variables
        y = {}
        n_ij = {}
        n_z = {}
        m={}
        a={}
        b={}
        C_st={}
        s={}
        for j in range(1,NA+1):
            a[j]=model.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
            b[j]=model.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
            C_st[j]=model.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
            s[j]=model.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
            for k in range(0,K+1):
                n_z[j,k]=model.addVar(vtype=GRB.BINARY, name = "n_z%s" % str([j,k]))
        for j in range(1,NA+1): 
            for i in range(1,I+1):
                y[i, j] = model.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
                n_ij[i,j]=model.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
        model.setObjective(quicksum(quicksum(k*n_z[j,k] for k in range(1,K+1)) for j in range (1,NA+1))
                           , GRB.MINIMIZE)

        for i in range(1,I+1):
            model.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
        for j in range(1,NA+1):
            model.addConstr(quicksum(n_z[j,k] for k in range(1,K+1))==1,"constriants 2")
        for j in range(1,NA+1):
            model.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
            for i in range(1,I+1):
                for k in range(0,K+1):
                    model.addConstr(n_ij[i,j]*k>=y[i,j]-Q[i]*(1-n_z[j,k]), "constriants 11")

        for j in range(1,NA):
            model.addConstr(quicksum(n_z[j,k]*k for k in range(1,K+1))<=
                            quicksum(n_z[j+1,k]*k for k in range(1,K+1)),"constriants 10") 
        model.addConstr(quicksum(quicksum(k*n_z[j,k] for k in range(1,K+1)) for j in range (1,NA+1))>=math.ceil(sum(Q.values())/n))
        model.addConstr(quicksum(quicksum(k*n_z[j,k] for k in range(1,K+1)) for j in range (1,NA+1))<=solutiondep[-1])
        if status[-1]!='none' and status[-1]!='infeasible' and NA<NAmax:
            startt=solution[-1]
            if startt[-1]-1<=K:
                startt[-1]=startt[-1]-1
                startt.append(1)
                startt=sorted(startt)
                
                for j in range(1,NA+1):
                    for k in range(1,K+1):
                        if k==startt[j-1]:
                            n_z[j,k].start=1
                #print (startt,NA,K,n_z)
        model.params.LogToConsole=0
        model.params.TimeLimit=120
        #model.write("zainali.lp")
        model.update()
        model.optimize()
        if model.status==GRB.INFEASIBLE:
            status.append('infeasible')
            solution.append('infeasible')
            solutiondep.append(solutiondep[-1])
        elif model.status!=GRB.OPTIMAL and model.MIPGap<1:
            solution3={j:sum(k*n_z[j,k].x for k in range(1,K+1)) for j in range(1,NA+1)}
            obj[NA]=Objective_value(list(solution3.values()))
            if sum(solution3.values())< solutiondep[-1]:
                solution.append(list(solution3.values()))
                solutiondep.append(sum(solution3.values()))
            else:
                solution.append(solution[-1])
                solutiondep.append(solutiondep[-1])
            status.append('unoptimal')
        elif model.status==GRB.OPTIMAL:
            status.append('optimal')
            solution3={j:sum(k*n_z[j,k].x for k in range(1,K+1)) for j in range(1,NA+1)}
            
            obj[NA]=Objective_value(list(solution3.values()))
            solution.append(list(solution3.values()))
            solutiondep.append(sum(solution3.values()))
            
            #if obj[NA]>obj[NA-1]:
                #solution.append(list(solution3.values()))
                #solutiondep.append(sum(solution3.values())) 
                #break
            #else:
            
            if max(solution3.values())==1:
                break
        else:
            
            obj[NA]=obj[NA-1]
            solution.append(solution[-1])
            solutiondep.append(solutiondep[-1])
            status.append('unoptimal')
        print (obj)
    return (solutiondep,obj,solution)


# In[ ]:





# In[3]:


####single deep rack model####
import math
def single_deep_model():
    #startTime=datetime.now()
    NA=int(math.ceil(sum(Q[i] for i in Q)/n))
    M = Model("MODEL")
    # Creating Variables
    y = {}
    n_ij = {}
    n_z = {}
    m={}
    a={}
    b={}
    C_st={}
    s={}
    for j in range(1,NA+1):
        a[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
        b[j]=M.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
        C_st[j]=M.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
    for j in range(1,NA+1): 
        for i in range(1,I+1):
            y[i, j] = M.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
    M.setObjective(quicksum(C_st[j] for j in range(1,NA+1))+C_v*NA+C_rw*w*h*l*n*NA, GRB.MINIMIZE)

    for i in range(1,I+1):
        M.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
    for j in range(1,NA+1):
        M.addConstr(quicksum(y[i,j] for i in range(1,I+1))<=n,"constriants 10")

    #M1
    for j in range(1,NA+1):
        M.addConstr(C_st[j]>=quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))*4/3*Txy*C_t,"constriants 9")
    M.params.LogToConsole=0
    #M.params.TimeLimit=60
    M.reset()
    M.update()
    M.optimize()
    return (M.ObjVal,NA)


# In[4]:


def model_3d(k):
    NA=NAmax
    m = Model("MODEL")
    y = {}
    n_ij = {}
    n_z = {}
    a={}
    b={}
    C_st={}
    s={}
    for j in range(1,NA+1):
        a[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
        b[j]=m.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
        C_st[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
        s[j]=m.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
        n_z[j]=m.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
    for j in range(1,NA+1): 
        for i in range(1,I+1):
            y[i, j] = m.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
            n_ij[i,j]=m.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
    m.setObjective(quicksum(C_st[j]+C_v*n_z[j]+C_rw*w*h*l*n*k*n_z[j] for j in range (1,NA+1))
                       ,GRB.MINIMIZE)
    m.addConstr(quicksum(n_z[j] for j in range(1,NA+1))>=sum(Q.values())/n/k ,"constriants 2")

    for i in range(1,I+1):
        m.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
    for j in range(1,NA):
        m.addConstr(n_z[j]>=n_z[j+1] ,"constriants 2")
    for j in range(1,NA+1):
        m.addConstr(sum(Q.values())*n_z[j]>=quicksum(y[i,j] for i in range (1,I+1)), 'cons')
        m.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
        for i in range(1,I+1):
            m.addConstr(n_ij[i,j]*k>=y[i,j], "constriants 11")
    for j in range(1,NA+1):

        m.addConstr(a[j]>=C_s+(1/n*EC31[k]+(1-1/n)*EC32[k])*C_t*
                    quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))-M3*(1-n_z[j]) ,"new constriants 3")
        m.addConstr(b[j]>=C_s*n+EC31[k]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                    -M3*(1-n_z[j]),"new constriants 4")
        m.addConstr(a[j]-b[j]<=M3*s[j],"new constriants 5")
        m.addConstr(b[j]-a[j]<=M3*(1-s[j]),"new constriants 6")
        m.addConstr(C_st[j]>=b[j]-M3*(1-s[j]),"constriants 7")
        m.addConstr(C_st[j]>=a[j]-M3*s[j],"constriants 8")
    m.params.LogToConsole=0
    m.params.TimeLimit=120
    m.update()
    m.optimize()
    
    if m.status==GRB.INFEASIBLE:
        return (math.inf,0,0)
    else:
        solution4={}
        for j in range(1,NA+1):
            if n_z[j].x==1:
                if s[j].x==0:
                    solution4[j]=1
                else:
                    solution4[j]=n
        return (m.ObjVal,k,sum(n_z[j].x for j in range(1,NA+1)))
def compact():
    objval=[]
    solution=[]
    for k in range(2,K+1):
        temp=model_3d(k)
        solution.append(temp)
        objval.append(temp[0])
        #print (solution) 
    print (solution)
    print (objval)
    return solution[objval.index(min(objval))]


# In[5]:


def rootnode_3d(NA):
    startTime=datetime.now()
    m = Model("MODEL")
    y = {}
    n_ij = {}
    n_z = {}
    K=max(Q.values())
    a={}
    b={}
    C_st={}
    s={}
    for j in range(1,NA+1):
        a[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "a%s" % str([j]))
        b[j]=m.addVar(lb=0,vtype=GRB.CONTINUOUS, name = "b%s" % str([j]))
        C_st[j]=m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = "C_st%s" % str([j]))
        s[j]=m.addVar(vtype=GRB.BINARY, name = "s%s" % str([j]))
        for k in range(0,K+1):
            n_z[j,k]=m.addVar(vtype=GRB.BINARY, name = "n_z%s" % str([j,k]))
    for j in range(1,NA+1): 
        for i in range(1,I+1):
            y[i, j] = m.addVar(lb=0, ub=Q[i],vtype=GRB.INTEGER, name = "y%s" % str([i,j]))
            n_ij[i,j]=m.addVar(lb=0, ub=n,vtype=GRB.INTEGER, name = "n_ij%s" % str([i,j]))
    m.setObjective(quicksum(C_st[j]+C_v*quicksum(n_z[j,k] for k in range(1,K+1))
                                +C_rw*w*h*l*n*quicksum(k*n_z[j,k] for k in range(0,K+1)) for j in range (1,NA+1))
                       , GRB.MINIMIZE)

    for i in range(1,I+1):
        m.addConstr(Q[i]==quicksum(y[i,j] for j in range (1,NA+1)),"constriants 1")
    for j in range(1,NA+1):
        m.addConstr(quicksum(n_z[j,k] for k in range(0,K+1))==1,"constriants 2")
    for j in range(1,NA+1):
        m.addConstr(quicksum(n_ij[i,j] for i in range(1,I+1))<=n,"constriants 10")
        for i in range(1,I+1):
            for k in range(0,K+1):
                m.addConstr(n_ij[i,j]*k>=y[i,j]-Q[i]*(1-n_z[j,k]), "constriants 11")
    for j in range(1,NA+1):
        for k in range(2,K+1):
            if k>=k1:
                m.addConstr(C_st[j]>=C_s*n+EC31[k]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                                -M3*(1-n_z[j,k]),"new constriants 1")
            elif k<=k2:
                m.addConstr(C_st[j]>=C_s+(1/n*EC31[k]+(1-1/n)*EC32[k])*C_t*
                            quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))-M3*(1-n_z[j,k]) ,"new constriants 2")
            else: 
                m.addConstr(a[j]>=C_s+(1/n*EC31[k]+(1-1/n)*EC32[k])*C_t*
                            quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))-M3*(1-n_z[j,k]) ,"new constriants 3")
                m.addConstr(b[j]>=C_s*n+EC31[k]*C_t*quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))
                            -M3*(1-n_z[j,k]),"new constriants 4")
        m.addConstr(a[j]-b[j]<=M3*s[j],"new constriants 5")
        m.addConstr(b[j]-a[j]<=M3*(1-s[j]),"new constriants 6")
        m.addConstr(C_st[j]>=b[j]-M3*(1-s[j]),"constriants 7")
        m.addConstr(C_st[j]>=a[j]-M3*s[j],"constriants 8")
        m.addConstr(C_st[j]>=quicksum(y[i,j]*p[i]/Q[i] for i in range(1,I+1))*4/3*Txy*C_t,"constriants 9")
    for j in range(1,NA):
        m.addConstr(quicksum(n_z[j,k]*k for k in range(0,K+1))<=
                        quicksum(n_z[j+1,k]*k for k in range(0,K+1)),"constriants 12")
        m.addConstr(n_z[j,1]==0,"constriants 12")

    #m.params.LogToConsole=0
    if NA>1:
        m.params.SolutionLimit=10
        m.params.TimeLimit=60
    #model.getEnv().set(GRB_IntParam_OutputFlag, 0)
   # model.write("zainali.lp")
    m.update()
    m.optimize()
    
    #print ((datetime.now() - startTime).total_seconds())
    solution1={(i,j):y[i,j].x for i in range(1,I+1) for j in range(1,NA+1)}
    solution2={(i,j):n_ij[i,j].x for i in range(1,I+1) for j in range(1,NA+1)}
    solution3={j:sum(k*n_z[j,k].x for k in range(0,K+1)) for j in range(1,NA+1) }
    solution4={j:0 for j in range(1,NA+1)}
    for j in range(1,NA+1):
        if s[j].x==0:
            for k in range(2,1+K):
                if n_z[j,k].x==1:
                    solution4[j]=1
        else:
            for k in range(2,1+K):
                if n_z[j,k].x==1:
                    solution4[j]=n
    print (solution3)
    return (m.ObjVal,solution1,solution2,solution3,solution4)


# In[14]:


from datetime import datetime # To calculate solution time
from gurobipy import *
import random as r
import math
import numpy as np
import itertools
import os
import xlwt
import xlrd 
import csv 
import pandas as pd# import this moduel so that i can read data from a file.
import ast
w=1.4
h=1.4
l=2
V_s=1.5
V_cx=2.5
V_cy=0.5
           
                                                                                                                                                                                                                                                                                                                                                                                
#######
NAmax=10
NAmin=1
path = '/Users/wenquandong/OneDrive - University of Tennessee/ASRS/system comparision and selection/experiment 5'
MYSGS = 'B&B with model mixed_p&Q new cost 125-2022'
newfolder = os.path.join(path, MYSGS)
#print (newfolder)
#os.makedirs(newfolder)
C_s=4
C_v=130
C_t=125
C_rw=0.1
I=27
Ni=0.2
Nd=0.0
Nd=0
savepath=os.path.join(newfolder,'model compare.xls'  )
wb=xlwt.Workbook()

readpath1=os.path.join(newfolder,'simulation data_time parameters.xls' )
readpath2=os.path.join(newfolder,'algorithm compare.xls' )
ps = xlrd.open_workbook(readpath1)
rs = xlrd.open_workbook(readpath2)
for f in range(4):
    columns=['obj_multi','solution_multi','obj_2d', 'solution_2d','obj_3d', 'solution_3d','shuttle_3d']
    Nd+=0.2
    summary=wb.add_sheet('summary %s'%(Nd))
    for i in range(7):                                                               
        summary.write(0,i,columns[i])
    result = rs.sheet_by_name('compare %s'%(Nd))
    for nq in range(5):
        solution = ast.literal_eval(result.cell_value(nq+1,3))
        parameter = ps.sheet_by_name('parameter%s%s'%(f,nq))
        q = parameter.col_values(1,1)
        P=parameter.col_values(2,1)
        #rint(P)
        I=len(q)
        n=30
        C_t=parameter.cell_value(1,5)
        Q={i:int(q[i-1]) for i in range(1,I+1)}
        p={i:float(P[i-1]) for i in range(1,I+1)}
        Txy=math.sqrt(w*h*n/V_cx/V_cy)
        K=max(Q.values())
        EC31={} #EC3Dj part 1 without m/n
        EC32={} #EC3Dj part 2
        EC31[0]=0
        EC32[0]=0
        EC31[1]=4/3*Txy
        EC32[1]=4/3*Txy
        for k in range(2,sum(Q.values())+1):
            EC32[k]=9/5*Txy+k*l/V_s
            if Txy>=2*k*l/V_s:
                EC31[k]=4/3*Txy+(2*k*l/V_s/Txy)**3/12*Txy
            else:
                EC31[k]=2/3*Txy+2*k*l/V_s*(1/2+(Txy/(2*k*l/V_s))**2/4)
        M1=sum(Q[i] for i in range(1,I+1))
        M2=max(w*h*sum(Q.values())/V_cx/V_cy,2*K*l/V_s)
        M3=C_s*n+(1/n*EC31[max(Q.values())]+(1-1/n)*EC32[max(Q.values())])*C_t 
        S=Objective_value(solution)
        #print (S)
        for j in range(len(S)):
            summary.write(nq+1,j,str(S[j]))
        single=single_deep_model()
        s_3d=compact()
        summary.write(nq+1,len(S),str(single[0]))
        summary.write(nq+1,len(S)+1,str(single[1]))
        summary.write(nq+1,len(S)+2,str(s_3d[0]))
        summary.write(nq+1,len(S)+3,str(s_3d[1]))
        summary.write(nq+1,len(S)+4,str(s_3d[2]))
        
wb.save(savepath)
   

