//root includes
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include <cmath>

using namespace std;

 void compute_stvs( const TVector3& p3mu, const TVector3& p3p, float& delta_pT, float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn )

{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y()) / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()- p3mu.Y()*delta_pT_vec.Y()) / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );
} 

void stv()
{
    //function to compute STVs from analyzer.C
   
//create the canvas
TCanvas	*c1= new TCanvas;
c1-> SetBottomMargin(0.5); //need some number, comeback
c1-> SetLeftMargin(0.5);

TH1D* f = new TH1D("STV_dist","STV Distribution", 50, -0.5, 15 ) ;

//import data from gst file
TFile* imp = new TFile("MicroBooNE_G18_10a_02_11a_numu_CCinclMEC.gst.root","read");
TTree* my_ptr= 0;
imp->GetObject("gst",my_ptr); //assign pointer of the tree

//set address for pointer of muon's momentum
double pxl =0, pyl =0, pzl = 0; //initialize the address
my_ptr->SetBranchAddress("pxl", &pxl);
my_ptr->SetBranchAddress("pyl", &pyl);
my_ptr->SetBranchAddress("pzl", &pzl);

long numb_ent = my_ptr -> GetEntries(); //get the number of entries

TVector3 p3mu; //initilize 3D vector of muon momentum
/*for (int i=0; i<= numb_ent; i++)
{
    my_ptr->GetEntry(i);
    p3mu.SetX(pxl);
    p3mu.SetY(pyl);
    p3mu.SetZ(pzl);
    f-> Fill(p3mu);
}*/
//f->Draw();

//begin of part 2
TVector3 p3p; //initilie 3D vector of leading proton momentum.
double pxf[100] , pyf[100], pzf[100] ; //array of momentum of particles in hadronic system.
my_ptr->SetBranchAddress("pxf",&pxf);
my_ptr->SetBranchAddress("pyf",&pyf);
my_ptr->SetBranchAddress("pzf",&pzf);

int nf =0; //number of final state particles in the hadronic system
my_ptr->SetBranchAddress("nf",&nf);
int pdgf[1000]; //PDG code of particle in hadronic system
my_ptr->SetBranchAddress("pdgf",&pdgf);
//loop to get value for pxl, pyl and pzl to compute p3mu
for (int i=0; i<= numb_ent; i++)
{
    my_ptr->GetEntry(i);
    p3mu.SetX(pxl);
    p3mu.SetY(pyl);
    p3mu.SetZ(pzl);
    bool isProton= false;
    double mag_p =0 ;// magnitude of proton momentum 
    int a =0; //index to keep track of highest momentum proton 
    for(int j=0; j<= nf; j++)
    {   
        
        if (pdgf[j] == 2212)
        {
            isProton= true;
            p3p.SetX(pxf[j]);
            p3p.SetY(pyf[j]);
            p3p.SetZ(pzf[j]);
            
            if (p3p.Mag()>mag_p) //only take the bigger value
            {
                mag_p = p3p.Mag(); 
                a = j;  //a is the index of the event which has the highest momentum 
                           
            };   
            
        };
    }
    p3p.SetX(pxf[a]);
    p3p.SetY(pyf[a]);
    p3p.SetZ(pzf[a]); 
    
    float delta_pT, delta_phiT, delta_alphaT, delta_pL, pn;
    compute_stvs(p3mu,p3p,delta_pT,delta_phiT,delta_alphaT,delta_pL, pn );
    
    if (isProton) f->Fill(delta_pT);

    f->Draw(); 




}







}






