#include "functions.h"


// containers
vector<graph> graphSeq;
vector<string> graphs;


// functions
void ReadIn(string name){
    double st = get_wall_time();
    cout << "\nRead in..." << endl;

    string fileName = "./graphSequences/" + name + "/" + "info.txt";

    string temp;
    ifstream fin(fileName);

    while(fin >> temp){
        graphs.push_back(temp);
    }

    for(int i = 0; i < graphs.size(); i++) {
        string fileName = "./graphSequences/" + name + "/" + graphs[i];
        graph pGraph(fileName);
        graphSeq.push_back(pGraph);
    }

    cout << "Instance: " << name << endl;
    cout << "Length of sequence: " << graphs.size() << endl;
    cout << "Processing time: " <<  fixed << setprecision(2) << (get_wall_time() - st) << " sec\n";
}


void GSIP_F2(int tau){
    double st = get_wall_time();
    int T = (int)graphSeq.size();
    cout << "\nGSIP-F2..." << endl;
    cout << "Find " << tau << "-persistent clique signature " << "T = " << graphs.size() << ", tau = " << tau << "..." << endl;

    GRBEnv* env = 0;
    GRBVar* xvar = 0;
    GRBVar* yvar = 0;
    GRBVar* zvar = 0;

    try{
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        //create variables
        xvar = model.addVars(graphSeq[0].n, GRB_BINARY);
        yvar=model.addVars((int)(graphSeq.size() - tau + 1), GRB_BINARY); // this variable is introduced to linearize GSIP-F2
        zvar=model.addVars((int)graphSeq.size(), GRB_BINARY);
        model.update();

        //create objective function, set bounds, name variables
        for(int i = 0; i < graphSeq[0].n; i++){
            xvar[i].set(GRB_DoubleAttr_Obj,1);
        }
        model.update();

        //create constraints
        for(int t = 0; t < (int)graphSeq.size(); t++){
            for(int i = 0; i < graphSeq[t].n; i++){
                for(int j = i + 1; j < graphSeq[t].n; j++){
                    if(graphSeq[t].IsAdj(i,j) == 0){
                        model.addConstr(xvar[i] + xvar[j] <= 2 - zvar[t]);
                    }
                }
            }
        }

        GRBLinExpr sumY = 0;
        for(int t = 0; t < (int)(graphSeq.size() - tau + 1); t++){
            sumY += yvar[t];
            for(int j = t; j <= t + tau - 1 ; j++){
                model.addConstr(yvar[t] <= zvar[j]);
            }
        }
        model.addConstr(sumY >= 1);
        model.update();

        // set Gurobi Parameters

        //set feasibility vs optimality balance
        model.set(GRB_IntParam_MIPFocus,0);
        //1-feasible sols quickly;2-prove optimality;3-focus on MIP bound; default is 0

        //set threads; review guidance on Gurobi.com; 0-default;
        model.set(GRB_IntParam_Threads,1);

        //set root node LPR solver
        model.set(GRB_IntParam_Method,-1);
        //-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent

        //set BC node LPR solver
        model.set(GRB_IntParam_NodeMethod,1);
        //0=primal simplex, 1=dual simplex, 2=barrier

        //set global cut aggressiveness; over-ridden by individual cut settings
        model.set(GRB_IntParam_Cuts,-1);
        //0=no cuts;1=moderate;2=aggressive;3=very aggressive;-1=default
        //model.set(GRB_IntParam_RLTCuts, 0);

        //set maximum time limit
        model.set(GRB_DoubleParam_TimeLimit,7200);

        //set termination gap limit; as needed; default is 1e-4
        model.set(GRB_DoubleParam_MIPGap,0);

        //set Gurobi log file name, if necessary; "" to switch off
        model.set(GRB_StringParam_LogFile, "");

        //set Gurobi screen display flag
        model.set(GRB_IntParam_OutputFlag,1);
        //0=switch off; 1=default

        // set Model Attributes

        //set objective to maximize
        model.set(GRB_IntAttr_ModelSense,-1);

        //set model name
        model.set(GRB_StringAttr_ModelName, "GSIP-F2");

        //in case of exhausting memory
        //model.getEnv().set(GRB_DoubleParam_NodefileStart,0.1);

        //begin optimization
        model.optimize();

        //output results
        if((int)model.get(GRB_IntAttr_SolCount) == 0){
            cout << "\nNo solution found, Gurobi optimization status: " << model.get(GRB_IntAttr_Status) << endl;

//            fstream fout("./graphSequences/" + seqName + "/results.txt", fstream::out);
//            fout << "No solution found, Gurobi optimization status: " << model.get(GRB_IntAttr_Status) << endl;
//            fout << "Processing time: " <<  fixed << setprecision(2) << (get_wall_time() - st) << " sec" << endl;
//            fout.close();
        }else{
            cout << "\nBest clique signature found: ";
            int count = 0;
            for(int i = 0; i < graphSeq[0].n; i++){
                if(xvar[i].get(GRB_DoubleAttr_X) > 0.5){
                    cout << i + 1 << " ";
                    count++;
                }
            }
            cout << endl;
            cout << "Clique signature size: " << count << endl;
            cout << "Best window: ";
            if(count != 0){
                for(int t = 0; t < graphSeq.size(); t++){
                    if(zvar[t].get(GRB_DoubleAttr_X) > 0.5){
                        cout << t + 1 << " (";
                        for(int i = t; i < t + tau; i++){
                            cout << "G" << i + 1;
                            if(i != t + tau -1){
                                cout << " ";
                            }
                        }
                        cout << ")" << endl;
                        break;
                    }
                }
            }else{
                cout << endl;
            }

            cout << "MIP gap: " << model.get(GRB_DoubleAttr_MIPGap) << endl;
            cout << "Best bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;

//            fstream fout("./graphSequences/" + seqName + "/results.txt", fstream::out);
//            fout << "Best clique signature found: ";
//            for(int i = 0; i < graphSeq[0].n; i++){
//                if(xvar[i].get(GRB_DoubleAttr_X) > 0.5){
//                    fout << i + 1 << " ";
//                }
//            }
//            fout << endl;
//            fout << "Windows: ";
//            for(int t = 0; t < graphSeq.size(); t++){
//                if(zvar[t].get(GRB_DoubleAttr_X) > 0.5){
//                    fout << t + 1 << " ";
//                }
//            }
//            fout << endl;
//            fout << "Clique signature size: " << count << endl;
//            fout << "MIP gap: " << model.get(GRB_DoubleAttr_MIPGap) << endl;
//            fout << "Best bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
//            fout << "Processing time: " << "Processing time: " <<  fixed << setprecision(2) << (get_wall_time() - st) << " sec" << endl;
//            fout.close();
        }

    }catch (GRBException e){
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e .getMessage() << endl;
    }catch(...){
        cout << "Exception during optimization" << endl;
    }
    cout << "Processing time: " <<  fixed << setprecision(2) << (get_wall_time() - st) << " sec" << endl;
}


void MW(int tau){
    double st = get_wall_time();
    int T = (int)graphSeq.size();
    cout << "\nMW..." << endl;
    cout << "Find " << tau << "-persistent clique signature " << "T = " << graphs.size() << ", tau = " << tau << "..." << endl;
    vector<int> bestClique;
    vector<int> flagOptimal(T - tau + 1, 1);
    bestClique.clear();
    int peel = 0;
    int bestWindow = 0;

    for(int i = 0; i < T - tau + 1; i++){
        cout << "\nIn window " << i+1 << "..." << endl;
        vector<int> pClique = GetPersistentCliqueWindow(i, tau, &peel, &flagOptimal);
        if(bestClique.size() < pClique.size()){
            bestClique = pClique;
            bestWindow = i;
        }
        if(i > 0){
            if(flagOptimal[i] == 0){
                if(flagOptimal[i - 1] == 0){
                    break;
                }
            }
        }
    }

    cout << "\nBest persistent clique found: ";
    for(int i = 0; i < bestClique.size(); i++){
        cout << bestClique[i] + 1 << " ";
    }
    cout << endl;
    cout << "Clique size: " << bestClique.size() << endl;
    if((int)bestClique.size() == 0){
        cout << "Best window: " << endl;
    }else{
        cout << "Best window: " << bestWindow + 1 << " (";
        for(int i = bestWindow; i < bestWindow + tau; i++){
            cout << "G" << i + 1;
            if(i < bestWindow + tau - 1){
                cout << " ";
            }
        }
        cout << ")" << endl;
    }
    cout << "Processing time: " <<  fixed << setprecision(2) << (get_wall_time() - st) << " sec" << endl;
}


vector<int> GetPersistentCliqueWindow(int windowHead, int tau, int* peelP, vector<int>* flagOptimalP){
    vector<int> bestClique;
    int upperBound = graphSeq[windowHead].n + 1;

    //create an intersection graph
    graph pGraph = GetIntersectionGraph(windowHead, tau);

    //get a heuristic maximal clique as lower bound (only in the first window)
    vector<int> heuClique;
    if(windowHead == 0){
        heuClique = GetHeuristicClique(&pGraph);
        *peelP = (int)heuClique.size();
    }

    //core peel
    CorePeel(&pGraph, peelP);

    //community peel
    CommunityPeel(&pGraph, peelP);

    //get components of size at least size of *peelP
    vector<vector<int>> components = GetComponents(&pGraph, peelP);

    //get vertex set after preprocessing
    vector<int> nodes;
    for(int i = 0; i < components.size(); i++){
        for(int j = 0; j < components[i].size(); j++){
            nodes.push_back(components[i][j]);
        }
    }

    if(nodes.size() > heuClique.size()){
        //solve
        GRBEnv* env=0;
        GRBVar* xvar=0;
        GRBVar* yvar=0;

        try{
            env = new GRBEnv();
            GRBModel model = GRBModel(*env);

            xvar = model.addVars(nodes.size(), GRB_BINARY);
            yvar = model.addVars((int)components.size(), GRB_BINARY);

            for(int i = 0; i < nodes.size(); i++){
                xvar[i].set(GRB_DoubleAttr_Obj, 1);
            }
            model.update();

            GRBLinExpr sumY = 0;
            int k = 0;
            for(int i = 0; i < components.size(); i++){
                sumY += yvar[i];
                for(int j = 0; j < components[i].size(); j++){
                    model.addConstr(xvar[k] <= yvar[i]);
                    k++;
                }
            }
            model.addConstr(sumY == 1);
            for(int i = 0; i < components.size(); i++){
                int begin = 0, end = 0;
                for(int j = 0; j < i; j++){
                    begin += (int)components[j].size();
                }
                end = begin + (int)components[i].size();

                for(int u = begin; u < end; u++){
                    for(int v = u + 1; v < end; v++){
                        if(pGraph.IsAdj(nodes[u], nodes[v]) == 0){
                            model.addConstr(xvar[u] + xvar[v] <= 1);
                        }
                    }
                }
            }
            model.update();

            // set Gurobi Parameters

            //set feasibility vs optimality balance
            model.set(GRB_IntParam_MIPFocus,0);
            //1-feasible sols quickly;2-prove optimality;3-focus on MIP bound; default is 0

            //set threads; review guidance on Gurobi.com; 0-default;
            model.set(GRB_IntParam_Threads,0);

            //set root node LPR solver
            model.set(GRB_IntParam_Method,-1);
            //-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent

            //set BC node LPR solver
            model.set(GRB_IntParam_NodeMethod,1);
            //0=primal simplex, 1=dual simplex, 2=barrier

            //set global cut aggressiveness; over-ridden by individual cut settings
            model.set(GRB_IntParam_Cuts,0);
            //0=no cuts;1=moderate;2=aggressive;3=very aggressive;-1=default

            //set maximum time limit
            model.set(GRB_DoubleParam_TimeLimit,3600);

            //set termination gap limit; as needed; default is 1e-4
            model.set(GRB_DoubleParam_MIPGap,0);

            //set Gurobi log file name, if necessary; "" to switch off
            model.set(GRB_StringParam_LogFile, "");

            //set Gurobi screen display flag
            model.set(GRB_IntParam_OutputFlag,1);
            //0=switch off; 1=default

            // set Model Attributes

            //set objective to maximize
            model.set(GRB_IntAttr_ModelSense,-1);

            //set model name
            model.set(GRB_StringAttr_ModelName, "FindPersistentClique");

            //in case of exhausting memory
            //model.getEnv().set(GRB_DoubleParam_NodefileStart,0.1);

            //begin optimization
            model.optimize();

            //get results
            if((int)model.get(GRB_IntAttr_SolCount)){
                for(int i = 0; i < nodes.size(); i++){
                    if(xvar[i].get(GRB_DoubleAttr_X) > 0.5){
                        bestClique.push_back(nodes[i]);
                    }
                }
                upperBound = model.get(GRB_DoubleAttr_ObjBound);
            }

        }catch (GRBException e){
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }catch(...){
            cout << "Exception during optimization" << endl;
        }
    }else{
        upperBound = heuClique.size();
    }

    if(heuClique.size() > bestClique.size()){
        if(*peelP < heuClique.size()){
            *peelP = (int)heuClique.size();
        }
        if(upperBound - heuClique.size() >= 1){
            (*flagOptimalP)[windowHead] = 0;
        }
        return heuClique;
    }else{
        if(*peelP < bestClique.size()){
            *peelP = (int)bestClique.size();
        }
        if(upperBound - bestClique.size() >= 1){
            (*flagOptimalP)[windowHead] = 0;
        }
        return bestClique;
    }
}

graph GetIntersectionGraph(int windowHead, int tau){
    if(tau == 1){
        graph pGraph = graphSeq[windowHead];
        return pGraph;
    }else{
        graph pGraph;
        pGraph.n = graphSeq[windowHead].n;

        for(int i = 0; i < pGraph.n; i++){
            node pNode(i);
            pGraph.nodeList.push_back(pNode);
        }

        int edgeCount = 0;
        for(int i = 0; i < pGraph.n; i++){
            for(int j = 0; j < graphSeq[windowHead].nodeList[i].neighbors.size(); j++){
                int pNeighbor = graphSeq[windowHead].nodeList[i].neighbors[j];
                if(pNeighbor > i){
                    int graphCount = 1;
                    for(int t =  windowHead + 1; t <= windowHead + tau -1; t++){
                        if(graphSeq[t].IsAdj(i, pNeighbor)){
                            graphCount++;
                        }else{
                            graphCount = 0;
                            break;
                        }
                    }
                    if(graphCount > 0){
                        pGraph.nodeList[i].neighbors.push_back(pNeighbor);
                        pGraph.nodeList[pNeighbor].neighbors.push_back(i);
                        edgeCount++;
                    }
                }
            }
        }
        pGraph.m = edgeCount;

        int maxDeg = 0;
        int maxDegNode = -1;
        for(int i = 0; i < pGraph.n; i++){
            pGraph.nodeList[i].degree = (int)pGraph.nodeList[i].neighbors.size();
            if(pGraph.nodeList[i].neighbors.size() > maxDeg){
                maxDeg = (int)pGraph.nodeList[i].neighbors.size();
                maxDegNode = i;
            }
        }
        pGraph.maxDeg = maxDeg;
        pGraph.maxDegNode = maxDegNode;
        return pGraph;
    }
}


vector<int> GetHeuristicClique(graph* graphP){
    vector<int> heuClique;
    heuClique.push_back(graphP->maxDegNode);

    vector<int>* poolP = &graphP->nodeList[graphP->maxDegNode].neighbors;
    vector<int> tempDegrees(poolP->size(), 0);
    for(int i = 0; i < poolP->size(); i++){
        int iNode = (*poolP)[i];
        for(int j = i + 1; j < poolP->size(); j++) {
            int jNode = (*poolP)[j];
            if(graphP->IsAdj(iNode, jNode)){
                tempDegrees[i]++;
                tempDegrees[j]++;
            }
        }
    }

    int maxDeg, maxDegNode, maxDegNodePos;
    do{
        maxDeg = 0;
        maxDegNode = -1;
        maxDegNodePos = -1;
        for(int i = 0; i < poolP->size(); i++){
            if(tempDegrees[i] > maxDeg){
                maxDeg = tempDegrees[i];
                maxDegNode = (*poolP)[i];
                maxDegNodePos = i;
            }
        }

        if(maxDeg > 0){
            heuClique.push_back(maxDegNode);
            tempDegrees[maxDegNodePos] = 0;
            for(int i = 0; i < poolP->size(); i++){
                int pNode = (*poolP)[i];
                if(graphP->IsAdj(pNode, maxDegNode)){
                    if(tempDegrees[i] > 0){
                        tempDegrees[i]--;
                    }
                }else{
                    tempDegrees[i] = 0;
                }
            }
        }
    }while(maxDeg > 0);

    return heuClique;
}


void CorePeel(graph* graphP, int* peelP){
    if(graphP->maxDeg > 0){
        queue<int> Q;
        vector<int> F(graphP->n, 0);
        for(int i = 0; i < graphP->n; i++){
            if(graphP->nodeList[i].degree < *peelP){
                Q.push(i);
                F[i] = 1;
            }
        }

        while(Q.size()){
            int pNode = Q.front();
            Q.pop();
            graphP->nodeList[pNode].degree = 0;
            for(int i = 0; i < graphP->nodeList[pNode].neighbors.size(); i++){
                int pNeighbor = graphP->nodeList[pNode].neighbors[i];
                if(graphP->nodeList[pNeighbor].degree > 0){
                    graphP->nodeList[pNeighbor].degree--;
                    if(graphP->nodeList[pNeighbor].degree > 0){
                        if(graphP->nodeList[pNeighbor].degree < *peelP){
                            if(F[pNeighbor] == 0){
                                Q.push(pNeighbor);
                                F[pNeighbor] = 1;
                            }
                        }
                    }
                }
            }
        }

        vector<int> neighbors;
        int maxDeg = 0;
        int maxDegNode = -1;
        int edgeCount = 0;
        for(int i = 0; i < graphP->n; i++){
            neighbors.clear();
            if(graphP->nodeList[i].degree > 0){
                for(int j = 0; j < graphP->nodeList[i].neighbors.size(); j++){
                    int pNeighbor = graphP->nodeList[i].neighbors[j];
                    if(graphP->nodeList[pNeighbor].degree > 0){
                        neighbors.push_back(pNeighbor);
                    }
                }
                graphP->nodeList[i].neighbors = neighbors;
                edgeCount += (int)neighbors.size();
                if(graphP->nodeList[i].degree > maxDeg){
                    maxDeg = graphP->nodeList[i].degree;
                    maxDegNode = i;
                }

            }else{
                graphP->nodeList[i].neighbors.clear();
            }
        }
        graphP->m = edgeCount/2;
        if(maxDeg > 0){
            graphP->maxDeg = maxDeg;
            graphP->maxDegNode = maxDegNode;
        }else{
            graphP->maxDeg = 0;
            graphP->maxDegNode = 0;
        }
    }
}


void CommunityPeel(graph* graphP, int* peelP){
    if(graphP->maxDeg > 0){
        queue<vector<int>> Q;
        for(int i = 0; i < graphP->n; i++){
            graphP->nodeList[i].adjacency.resize(graphP->n, 0);
            graphP->nodeList[i].F.resize(graphP->n, 0); // used as flag to indicate if an edge is put into Q or not
            graphP->nodeList[i].numCommonNeighbors.resize(graphP->n, 0);
            graphP->nodeList[i].commonNeighbors.resize(graphP->n);
        }

        for(int i = 0; i < graphP->n; i++){
            for(int j = 0; j < graphP->nodeList[i].neighbors.size(); j++){
                int pNeighbor = graphP->nodeList[i].neighbors[j];
                if(graphP->nodeList[i].adjacency[pNeighbor] == 0){

                    graphP->nodeList[i].adjacency[pNeighbor] = 1;
                    graphP->nodeList[pNeighbor].adjacency[i] = 1;

                    if(i < pNeighbor){
                        graphP->nodeList[i].commonNeighbors[pNeighbor] = graphP->FindCommonNeighbors(i, pNeighbor);
                        graphP->nodeList[i].numCommonNeighbors[pNeighbor] = (int)graphP->nodeList[i].commonNeighbors[pNeighbor].size();
                        graphP->nodeList[pNeighbor].numCommonNeighbors[i] = graphP->nodeList[i].numCommonNeighbors[pNeighbor];
                    }else{
                        graphP->nodeList[pNeighbor].commonNeighbors[i] = graphP->FindCommonNeighbors(i, pNeighbor);
                        graphP->nodeList[pNeighbor].numCommonNeighbors[i] = (int)graphP->nodeList[pNeighbor].commonNeighbors[i].size();
                        graphP->nodeList[i].numCommonNeighbors[pNeighbor] = graphP->nodeList[pNeighbor].numCommonNeighbors[i];
                    }

                    if(graphP->nodeList[i].numCommonNeighbors[pNeighbor] < ((*peelP) - 1)){
                        if(graphP->nodeList[i].F[pNeighbor] == 0){
                            vector<int> tempVec{i, pNeighbor};
                            Q.push(tempVec);
                            graphP->nodeList[i].F[pNeighbor] = 1;
                            graphP->nodeList[pNeighbor].F[i] = 1;
                        }
                    }
                }
            }
        }

        while(Q.size()){
            vector<int> pEdge = Q.front();
            Q.pop();
            int u = pEdge[0], v = pEdge[1];
            int i, j;
            if(u < v){
                i = u;
                j = v;
            }else{
                i = v;
                j = u;
            }

            graphP->nodeList[i].adjacency[j] = 0;
            graphP->nodeList[j].adjacency[i] = 0;
            for(int l = 0; l < graphP->nodeList[i].commonNeighbors[j].size(); l++){
                int pCommonNeighbor = graphP->nodeList[i].commonNeighbors[j][l];

                if(graphP->nodeList[i].adjacency[pCommonNeighbor] > 0 && graphP->nodeList[j].adjacency[pCommonNeighbor] > 0){
                    graphP->nodeList[i].numCommonNeighbors[pCommonNeighbor]--;
                    graphP->nodeList[pCommonNeighbor].numCommonNeighbors[i]--;
                    if(graphP->nodeList[i].numCommonNeighbors[pCommonNeighbor] < ((*peelP) - 1)){
                        if(graphP->nodeList[i].F[pCommonNeighbor] == 0){
                            vector<int> tempVec{i, pCommonNeighbor};
                            Q.push(tempVec);
                            graphP->nodeList[i].F[pCommonNeighbor] = 1;
                            graphP->nodeList[pCommonNeighbor].F[i] = 1;
                        }
                    }

                    graphP->nodeList[j].numCommonNeighbors[pCommonNeighbor]--;
                    graphP->nodeList[pCommonNeighbor].numCommonNeighbors[j]--;
                    if(graphP->nodeList[j].numCommonNeighbors[pCommonNeighbor] < ((*peelP) - 1)){
                        if(graphP->nodeList[j].F[pCommonNeighbor] == 0){
                            vector<int> tempVec{j, pCommonNeighbor};
                            Q.push(tempVec);
                            graphP->nodeList[j].F[pCommonNeighbor] = 1;
                            graphP->nodeList[pCommonNeighbor].F[j] = 1;
                        }
                    }
                }

            }
        }

        vector<int> neighbors;
        int maxDeg = 0;
        int maxDegNode = -1;
        int edgeCount = 0;
        for(int i = 0; i < graphP->n; i++){
            neighbors.clear();
            for(int j = 0; j < graphP->nodeList[i].neighbors.size(); j++){
                int pNeighbor = graphP->nodeList[i].neighbors[j];
                if(graphP->nodeList[i].adjacency[pNeighbor] == 1){
                    neighbors.push_back(pNeighbor);
                }
            }
            graphP->nodeList[i].neighbors = neighbors;
            graphP->nodeList[i].degree = (int)neighbors.size();
            edgeCount += graphP->nodeList[i].degree;
            if(maxDeg < graphP->nodeList[i].degree){
                maxDeg = graphP->nodeList[i].degree;
                maxDegNode = i;
            }
        }
        graphP->m = edgeCount/2;
        if(maxDeg > 0){
            graphP->maxDeg = maxDeg;
            graphP->maxDegNode = maxDegNode;
        }else{
            graphP->maxDeg = 0;
            graphP->maxDegNode = 0;
        }

        for(int i = 0; i < graphP->n; i++){
            vector<int>().swap(graphP->nodeList[i].adjacency);
            vector<int>().swap(graphP->nodeList[i].F);
            vector<int>().swap(graphP->nodeList[i].numCommonNeighbors);
            vector<vector<int>>().swap(graphP->nodeList[i].commonNeighbors);
        }
    }
}


vector<vector<int>> GetComponents(graph* graphP, int* peelP){
    vector<vector<int>> components;
    if(graphP->maxDeg >= *peelP){
        vector<int> F(graphP->n, 0);
        vector<int> component;
        for(int i = 0; i < graphP->n; i++){
            if(F[i] == 0){
                component.clear();
                F[i] = 1;
                component = graphP->BFS(i);
                component.push_back(i);
                sort(component.begin(), component.end());
                for(int j = 0; j < component.size(); j++){
                    int pNode = component[j];
                    F[pNode] = 1;
                }
                if((int)component.size() >= *peelP){
                    components.push_back(component);
                }
            }
        }

        vector<int>().swap(F);
        vector<int>().swap(component);
        sort(components.begin(), components.end(), cmp);//sort components in descending order of component size
    }
    return components;
}


//sort a vector of vector in descending order of vector size
bool cmp(const vector<int> &a,const vector<int> &b)
{
    return a.size() > b.size();
}


string itos_c(int i){
    stringstream s;
    s << i;
    return s.str();
}


int stoi_c(string i){
    stringstream geek(i);
    int s=0;
    geek >> s;
    return s;
}


double stod_c(string i){
    stringstream geek(i);
    double s=0;
    geek >> s;
    return s;
}


string dtos_c(double i){
    stringstream s;
    s << i;
    return s.str();
}


void makeDir(string dirName){
    mkdir(dirName.c_str(), S_IRWXU);
}


void emptyDir(string dirName){
    system(("rm -r " + dirName +"/*").c_str());
}


void goToDir(string dirName){
    chdir(dirName.c_str());
}


string getDir() {
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX);
    string current_working_dir(buff);
    return current_working_dir;
}


double get_wall_time(){
    struct timeval time;
    if(gettimeofday(&time, NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec*.000001;
}


double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


