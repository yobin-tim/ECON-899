#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include <quadpack.h>
#import <maximize>
#include <oxprob.h>
/**********************************************************************************************/
/* Global variables */
decl alpha=2;
decl lambda=-4;
decl beta=0.99;
decl ibar=8;
decl ps=1;
decl pr=4;
decl c=1;
decl gamma=1/2;
decl mF1,mF0,mS;
/* Simulated choices and states */
decl vPhat; /* Estimted choice-probalities Sx1 */
decl vY; /* Observed choices Nx1 */
decl vSid; /* Observed state IDs Nx1 */
/* State space indices */
enum{iI,iC,iP}; /* Columns identifiers for each state variable */
/* Parameter names */
enum{ilambda}; /* Row identifiers for parameter */

main(){
  format(1000);
  decl mPi=<0.75,0.25;0.9,0.1>;
  decl mGamma=(gamma~(1-gamma))|(gamma~(1-gamma));
  println("Gamma: ",mGamma);
  decl vIgrid=range(0,ibar,c)';
  decl vPgrid=pr|ps;
  decl vCgrid=0|c;
  decl S=rows(vIgrid)*rows(vCgrid)*rows(vPgrid);
  decl mSid=(vIgrid|vIgrid|vIgrid|vIgrid)
    ~(zeros(rows(vIgrid),1)|ones(rows(vIgrid),1)|zeros(rows(vIgrid),1)|ones(rows(vIgrid),1))
    ~(zeros(2*rows(vIgrid),1)|ones(2*rows(vIgrid),1));
  mS=mSid;
  mS[][iC]=vCgrid[mS[][iC]];
  mS[][iP]=vPgrid[mS[][iP]];
  println("State space: ","%c",{"I","C","P"},mS);
  savemat("PS4_state_space.csv",range(0,S-1)'~mS,{"id","I","C","P"});
  println("State points: ",S);

  decl vU1=alpha*mS[][iC]-mS[][iP];
  decl vU0=alpha*mS[][iC].*(mS[][iI].>0)+lambda*(mS[][iC].>0).*(mS[][iI].==0);
  println("%c",{"I","C","P","U1","U0"},mS~vU1~vU0);

  decl i,j;
  mF0=new matrix[S][S];
  mF1=new matrix[S][S];
  decl vI1,vI0;
  vI1=setbounds(mS[][iI]-mS[][iC]+1,0,ibar);
  vI0=setbounds(mS[][iI]-mS[][iC],0,ibar);
  decl aStateID=new array[S];
  for(j=0;j<S;j++)
    {
      mF0[][j]=mPi[mSid[][iP]][mSid[j][iP]].*mGamma[mSid[][iC]][mSid[j][iC]].*(mS[j][iI].==vI0);
      mF1[][j]=mPi[mSid[][iP]][mSid[j][iP]].*mGamma[mSid[][iC]][mSid[j][iC]].*(mS[j][iI].==vI1);
      aStateID[j]=sprint("S",j);
    }
  savemat("PS4_transition_a0.csv",range(0,S-1)'~mF0,{"id"}~aStateID);
  savemat("PS4_transition_a1.csv",range(0,S-1)'~mF1,{"id"}~aStateID);  

  
  /* CCP iteration */
  decl vP=constant(1/2,S,1),vP0,vEV;
  println(vP0);
}
