/*********************************************
 * OPL 12.8.0.0 Model
 * Author: ofers
 * Creation Date: Nov 16, 2022 at 14:02:19
 *********************************************/
int NumOfObj=...;
range rObjectives=1..NumOfObj;
int MaxNumOfPoints=...;
int n=...;
range Dimension=1..n;
float w[rObjectives]=...;
int LB=...;
int UB=...;
dvar float x[d in Dimension] in LB..UB;
dvar int z[d in Dimension] in LB..UB; // <=== INT
dvar float obj[r in rObjectives];
float epsilon=...;
dvar float alpha;
float walpha = 100000;
float C=...;
float solution[0..MaxNumOfPoints][rObjectives]=...;
int iteration=...;
int run_time = ...;
tuple betaType{
	int ro;
	int point;
}
{betaType} betaSet={<i,j>|i in rObjectives,j in 0..MaxNumOfPoints:j<=iteration};
dvar boolean beta[betaSet];
int L=2000; // <===== LARGE INTEGER
float c1[1..2*n]=...;
float c2[1..2*n]=...;
float Hess[1..2*n][1..2*n]=...;
execute {
    cplex.tilim = run_time; //computation time limit
    cplex.epagap = epsilon;
	cplex.epgap = epsilon;
	cplex.epint = epsilon;
	cplex.epopt = epsilon;
	cplex.eprhs = epsilon;
	cplex.polishafterintsol = 1;
	cplex.threads=8;
	cplex.mipdisplay=4;
}   
/******************************************************************************
 * MODEL
 ******************************************************************************/
dexpr float y[i in 1..2*n] = (i in Dimension)?x[i]:z[i-n];
minimize walpha*alpha + sum(m in rObjectives)(w[m]*obj[m]);
subject to{
	obj[1]==(sum(d1,d2 in 1..2*n) Hess[d1][d2]*(y[d1]-c1[d1])*(y[d2]-c1[d2]))/C; 
	obj[2]==(sum(d1,d2 in 1..2*n) Hess[d1][d2]*(y[d1]-c2[d1])*(y[d2]-c2[d2]))/C; 
	
	forall(m in rObjectives,p in 0..MaxNumOfPoints:p<=iteration){
		alpha>=obj[m]-solution[p][m]+(1-beta[<m,p>])*L;
	}
	forall(p in 0..MaxNumOfPoints:p<=iteration){
		sum(m in rObjectives)(beta[<m,p>])==1;
	}
}

float temp[i in 0..MaxNumOfPoints][j in rObjectives]=solution[i][j];
execute{
  var flag=0;
  for (var i in iteration) {
 	   	 if((temp[iteration+1][1]==obj[1]&&temp[iteration+1][2]==obj[2]&& temp[iteration+1][3]==obj[3])) { flag=1; }
  }
  if(flag==0) {
	for(var i in rObjectives){ temp[iteration+1][i]=obj[i]; }
	var ofile=new IloOplOutputFile("EPF.dat");
	var itr=iteration+1;
	ofile.writeln("iteration="+itr+";");
	ofile.writeln("solution="+temp+";");
	ofile.close();
  }	
}