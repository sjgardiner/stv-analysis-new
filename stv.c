#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include <cmath>
#include "TCanvas.h"
using namespace std;

//function to calculate STVs
void compute_stvs( const TVector3& p3mu, const TVector3& p3p, float& delta_pT, float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn )

{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y()) / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()- p3mu.Y()*delta_pT_vec.Y()) / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );
}



void stv2(){
//Canvas for muon momentum
TCanvas	*c1= new TCanvas;
TH1D* mu_hist = new TH1D("STV_dist","Muon Momentum- Default", 200, -0.5, 5 ) ;
TH1D* mu_hist2 = new TH1D("STV_dist","Muon Momentum- Alternate ", 200, -0.5, 5 ) ;

//import data from gst output file
TFile* imp = new TFile("MicroBooNE_G18_10a_02_11a_numu_CCinclMEC.gst.root","read"); //default model
TFile* imp2 = new TFile("MicroBooNE_G00_00b_00_000_numu_CCinclMEC.gst.root","read"); //alternate model
TTree* my_ptr= 0;
TTree* my_ptr2= 0;
imp->GetObject("gst",my_ptr); //point to gst_default
imp2->GetObject("gst",my_ptr2); // point to gst alternate

//set address for pointer of muon's momentum
double pxl = 0, pyl = 0, pzl = 0; //initialize the address
//default model
my_ptr->SetBranchAddress("pxl", &pxl);
my_ptr->SetBranchAddress("pyl", &pyl);
my_ptr->SetBranchAddress("pzl", &pzl);
double pxl2 = 0, pyl2 = 0, pzl2 = 0; //alternate model
my_ptr2->SetBranchAddress("pxl", &pxl2);
my_ptr2->SetBranchAddress("pyl", &pyl2);
my_ptr2->SetBranchAddress("pzl", &pzl2);
long numb_ent = my_ptr -> GetEntries(); //get the number of entries 10e6

TVector3 p3mu, p3mu2; //initilize 3D vector of muon momentum

double mpa, mpa2; //variable to store the value of muon polar angle
TH1D* mpa_hist = new TH1D("STV_dist","Muon polar angle- Default", 200, -1, 1 ) ; //TH1D for mpa default model
TH1D* mpa_hist2 = new TH1D("STV_dist","Muon polar angle- Alternate", 200, -1, 1 ) ; //TH1D for mps alternate model

int nf =0; //number of final state particles in the hadronic system
my_ptr->SetBranchAddress("nf",&nf);
int pdgf[1000]; //PDG code of particle in hadronic system
int pdgf2[1000]; //PDG code of particle in hadronic system for alternate model
my_ptr->SetBranchAddress("pdgf",pdgf);
my_ptr2->SetBranchAddress("pdgf",pdgf2);

TVector3 p3p, p3p2; //initilie 3D vector of leading proton momentum for both models

double pxf[100] , pyf[100], pzf[100] ; //array of momentum of particles in hadronic system.
my_ptr->SetBranchAddress("pxf",pxf);
my_ptr->SetBranchAddress("pyf",pyf);
my_ptr->SetBranchAddress("pzf",pzf);

double pxf2[100] , pyf2[100], pzf2[100] ; //array of momentum of particles in hadronic system for alternate model
my_ptr->SetBranchAddress("pxf",pxf2);
my_ptr->SetBranchAddress("pyf",pyf2);
my_ptr->SetBranchAddress("pzf",pzf2);

TH1D* leadp_hist = new TH1D("STV_dist","Leading proton momentum- Default", 200, -0.5, 5 ) ; //TH1D for leading proton -default
TH1D* leadp_hist2= new TH1D("STV_dist","Leading proton momentum- Alternate", 200, -0.5, 5 ) ;//TH1D for leading proton -alternate

double ppa, ppa2; //variable to store proton polar angle
TH1D* ppa_hist = new TH1D("STV_dist","proton polar angle- Default", 200, -1, 1 ) ; //TH1D for proton polar angle -default
TH1D* ppa_hist2= new TH1D("STV_dist","proton polar angle- Alternate", 200, -1, 1 ) ;//TH1D for proton polar angle -alternate

//float delta_pT, delta_pT2; //variable for delta pT in STVs
TH1D* pT_hist = new TH1D("STV_dist","delta pT- Default", 200, -0.5, 5 ) ;  //TH1D for delta pT-default
TH1D* pT_hist2= new TH1D("STV_dist","delta pT- Alternate", 200, -0.5, 5 ); //TH1D for delta pT -alternate

TH1D* alphaT_hist = new TH1D("STV_dist","delta alphaT- Default", 200, -0.5, 5 ) ;  //TH1D for delta alphaT-default
TH1D* alphaT_hist2= new TH1D("STV_dist","delta alphaT- Alternate", 200, -0.5, 5 ); //TH1D for delta alphaT-alternate

TH1D* phiT_hist = new TH1D("STV_dist","delta phiT- Default", 200, -0.5, 5 ) ;  //TH1D for delta phiT-default
TH1D* phiT_hist2= new TH1D("STV_dist","delta phiT- Alternate", 200, -0.5, 5 ); //TH1D for delta phiT-alternate


for (int i=0; i<= numb_ent; i++)
{
    my_ptr->GetEntry(i);   my_ptr2->GetEntry(i);

    p3mu.SetX(pxl);   p3mu2.SetX(pxl2);
    p3mu.SetY(pyl);   p3mu2.SetY(pyl2);
    p3mu.SetZ(pzl);   p3mu2.SetZ(pzl2);
    mu_hist-> Fill(p3mu.Mag());   //muon momentum calculation
    mu_hist2-> Fill(p3mu2.Mag());
//calculating muon polar angle
    mpa = pzl/ p3mu.Mag();  mpa2 = pzl2/ p3mu2.Mag();
    mpa_hist-> Fill(mpa);
    mpa_hist2-> Fill(mpa2);

//calculating momentum of leading proton
    bool isProton= false;
    double mag_p= 0, mag_p2= 0 ;// magnitude of proton momentum
    int a= 0, a2= 0; //index to keep track of highest momentum proton
    //loop to find highest momentum proton - default model
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
    p3p.SetX(pxf[a]); //set the coordinate of the highest momentum -default
    p3p.SetY(pyf[a]);
    p3p.SetZ(pzf[a]);
    cout << p3p.Mag() << endl;
    ppa = pzf[a] / p3p.Mag();  //calculate cosine of the proton polar angle
    float delta_pT, delta_phiT, delta_alphaT, delta_pL, pn;
    compute_stvs(p3mu,p3p,delta_pT,delta_phiT,delta_alphaT,delta_pL, pn );


   if(isProton)
    {
        leadp_hist->Fill(p3p.Mag()); //fill the magnitude of leading proton
        ppa_hist->Fill(ppa); //fill the cosine of proton polar angle
        pT_hist->Fill(delta_pT); //fill pT -default
        alphaT_hist->Fill(delta_alphaT); //fill alphaT-default
        phiT_hist->Fill(delta_phiT); //fill phiT-default
    }
    //for loop to find leading proton -alternate
    for(int k=0; k<= nf; k++)
    {

        if (pdgf2[k] == 2212)
        {
            isProton= true;
            p3p2.SetX(pxf[k]);
            p3p2.SetY(pyf[k]);
            p3p2.SetZ(pzf[k]);

            if (p3p2.Mag()>mag_p2) //only take the bigger value
            {
                mag_p2 = p3p.Mag();
                a2 = k;  //a is the index of the event which has the highest momentum

            };

        };

    }
    p3p2.SetX(pxf[a2]); // set the coordinate of the highest momentum -alternate
    p3p2.SetY(pyf[a2]);
    p3p2.SetZ(pzf[a2]);
    ppa2 = pzf[a2] / p3p2.Mag();  //calculate cosine of the proton polar angle
    float delta_pT2, delta_phiT2, delta_alphaT2, delta_pL2, pn2;
    compute_stvs(p3mu2,p3p2,delta_pT2,delta_phiT2,delta_alphaT2,delta_pL2, pn2 );

    if(isProton)
    {
        leadp_hist2->Fill(p3p2.Mag()); //fill the magnitude of leading proton
        ppa_hist2->Fill(ppa2); //fill the cosine of proton polar angle.
        pT_hist2->Fill(delta_pT2); //fill pT -alternate
        alphaT_hist2->Fill(delta_alphaT2); //fill alphaT-alternate
        phiT_hist2->Fill(delta_phiT2);// fill phiT-alternate
    }



}
mu_hist->Draw();
mu_hist2->Draw("SAME");

//Canvas for cosine the muon polar angle (mpa)
TCanvas *c2 =  new TCanvas;
mpa_hist->Draw();
mpa_hist2->Draw("SAME");

//Canvas for leading proton momentum
TCanvas *c3 =  new TCanvas;
leadp_hist->Draw();
leadp_hist2->Draw("SAME");

//Canvas for cosine of the proton polar angle
TCanvas *c4 =  new TCanvas;
ppa_hist->Draw();
ppa_hist2->Draw("SAME");

//Canvas for pT
TCanvas *c5 =  new TCanvas;
pT_hist->Draw();
pT_hist2->Draw("SAME");

//Canvas for alphaT
TCanvas *c6 =  new TCanvas;
alphaT_hist->Draw();
alphaT_hist2->Draw("SAME");

//Canvas for phiT
TCanvas *c7 =  new TCanvas;
phiT_hist->Draw();
phiT_hist2->Draw("SAME");















}
