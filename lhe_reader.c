// use this code in root
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// The above line will produce executable file name lhe_reader_non_decayed
// ./lhe_reader_non_decayed unweighted_events








#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "TCanvas.h"
using namespace std;


int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  
  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;
  //TH1::SetDefaultSumw2(true);
  TH1F** hPt_bq = new TH1F*[4]; 
  TH1F** hEta_bq = new TH1F*[4];
  TH1F** hPhi_bq = new TH1F*[4];

  TH1F** hPt_q = new TH1F*[4];
  TH1F** hEta_q = new TH1F*[4];
  TH1F** hPhi_q = new TH1F*[4];
   //cout<<"I am here: "<<endl;
  for (int iquark=0 ; iquark<4 ; iquark++) {
    ostringstream oss;
    oss << "q" << iquark ;
    string ptstr = oss.str()+"_pt";
    hPt_q[iquark] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_q[iquark] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_q[iquark] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.5,3.5) ;
  }

    for (int bquark=0 ; bquark<4 ; bquark++) {
    ostringstream oss;
    oss << "bq" << bquark ;
    string ptstr = oss.str()+"_pt";
    hPt_bq[bquark] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_bq[bquark] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_bq[bquark] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
 
  TH1F** hPt_mu = new TH1F*[4];
  TH1F** hEta_mu = new TH1F*[4];
  TH1F** hPhi_mu = new TH1F*[4];
  for (int mu=0 ; mu<4 ; mu++) {
    ostringstream oss;
    oss << "muon" << mu ;
    string ptstr = oss.str()+"_pt";
    hPt_mu[mu] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_mu[mu] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_mu[mu] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
  TH1F** hPt_elec = new TH1F*[4];
  TH1F** hEta_elec = new TH1F*[4];
  TH1F** hPhi_elec = new TH1F*[4];
  for (int ele=0 ; ele<4 ; ele++) {
    ostringstream oss;
    oss << "elec" << ele ;
    string ptstr = oss.str()+"_pt";
    hPt_elec[ele] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_elec[ele] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_elec[ele] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.16,3.16) ;
  }
  TH1F** hPt_wboson = new TH1F*[4];
  TH1F** hEta_wboson = new TH1F*[4];
  TH1F** hPhi_wboson = new TH1F*[4];

  for (int wboson=0 ; wboson<4 ; wboson++) {
    ostringstream oss;
    oss << "wboson" << wboson ;
    string ptstr = oss.str()+"_pt";
    hPt_wboson[wboson] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_wboson[wboson] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_wboson[wboson] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.16,3.16) ;
    string massstr=oss.str()+"_mass";
   
  }

  TH1F** hPt_top = new TH1F*[4];
  TH1F** hEta_top = new TH1F*[4];
  TH1F** hPhi_top = new TH1F*[4];
  TH1F** hM_top = new TH1F*[4];
  for (int t=0 ; t<4 ; t++) {
    ostringstream oss;
    oss << "top_quark" << t ;
    string ptstr = oss.str()+"_pt";
    hPt_top[t] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,600) ;
    string etastr=oss.str()+"_eta";
    hEta_top[t] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_top[t] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.5,3.5) ;
    string massstr=oss.str()+"_mass";
   
  }
int nlept=0, nsemi=0, nhadr=0;
  ifstream ff(lhefname.c_str(),ios::in); 
  int negativeWeight = 0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;


    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;

      event++;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }

      if (event%1==0) cout << "reading event "<< event << endl;
     
     
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
      int n_muon=0, n_elec=0, n_mu_nu=0, n_elec_nu=0, n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bbar=0;
      int  bj[2]={-1,-1};
      int muon[4]={-1,-1,-1,-1};
      int mu_nu[4]={-1,-1,-1,-1};
      int elec[4]={-1,-1,-1,-1};
      int elec_nu[4]={-1,-1,-1,-1};
      int bq[4]={-1,-1,-1,-1};
      int q[4]={-1,-1,-1,-1};
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_muon1, v_muon2, v_muon3, v_muon4,v_muon, v_elec, v_mu_nu1, v_mu_nu2, v_mu_nu3, v_mu_nu4, v_mu_nu, v_elec_nu;
      TLorentzVector v_elec1, v_elec2, v_elec3, v_elec4, v_elec_nu1, v_elec_nu2, v_elec_nu3, v_elec_nu4;
      TLorentzVector v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4;
      TLorentzVector v_bquark1, v_bquark2, v_bquark3, v_bquark4;
      TLorentzVector v_top1, v_top2, v_top3, v_top4;
      
       Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      spin[0]= -99999;
      for (int i=1 ; i<npart+1 ; i++) {
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        
        
        
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        
         if (abs(Id[Mother1[i]])==6) { 
	  if (abs(Id[i])==5) { 

	         if (bq[0] == -1) bq[0]=i;
            else if (bq[1]==-1) bq[1]=i;
            else if (bq[2]==-1) bq[2]=i;
            else if (bq[3]==-1) bq[3]=i;
            else cout << "ERROR : more than 4 b quarks" << endl;
           }
	}
        
	if (abs(Id[Mother1[i]])==24) { 
         if (abs(Id[i]) < 5) { 
            cout << q[i] << endl;
            if (q[0] == -1) q[0]=i;
            else if (q[1]==-1) q[1]=i;
            else if (q[2]==-1) q[2]=i;
            else if (q[3]==-1) q[3]=i;
            else cout << "ERROR : more than 4 light quarks" << endl;
          }

        if ( abs(Id[i])==13 ) { 
	         if (muon[0] == -1) muon[0]=i;
            else if (muon[1]==-1) muon[1]=i;
            else if (muon[2]==-1) muon[2]=i;
            else if  (muon[3]==-1) muon[3]=i;
	    else cout << "ERROR : more than 4 b muons" << endl;



          n_muon++;     
          }
          if ( abs(Id[i])==14 ) {//muon neutrinos
	         if (mu_nu[0] == -1) mu_nu[0]=i;
            else if (mu_nu[1]==-1) mu_nu[1]=i;
            else if (mu_nu[2]==-1) mu_nu[2]=i;
            else if (mu_nu[3]==-1) mu_nu[3]=i;
            else cout << "ERROR : more than 4 muon neutrinos" << endl;


          n_mu_nu++;
	  }

        if ( abs(Id[i])==11 ) { // electrons
	         if (elec[0] == -1) elec[0]=i;
            else if (elec[1]==-1) elec[1]=i;
            else if (elec[2]==-1) elec[2]=i;
            else if (elec[3]==-1) elec[3]=i;
	    else cout << "ERROR : more than 4 electrons" << endl;


          n_elec++;
          }
          if ( abs(Id[i])==12 ) {
	         if (elec_nu[0] == -1) elec_nu[0]=i;
            else if (elec_nu[1]==-1) elec_nu[1]=i;
            else if (elec_nu[2]==-1) elec_nu[2]=i;
            else if (elec_nu[3]==-1) elec_nu[3]=i;
	    else cout << "ERROR : more than 4 electron neutrinos" << endl;



          n_elec_nu++;
          }

	}

} 
  cout<<" no: of muon:  "<<n_muon<<",  no. of muon nu:  "<<n_mu_nu<<endl;
  cout<<" no: of elec:  "<<n_elec<<",  no. of elec nu:  "<<n_elec_nu<<endl;

double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};
	std::sort(bq_pt,bq_pt+4);
	double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
	std::sort(bq_eta,bq_eta+4);
	double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
	std::sort(bq_phi,bq_phi+4);

  	v_bquark1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
  	v_bquark2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
  	v_bquark3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
  	v_bquark4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
   
     for (int nbq=0; nbq<4 ; nbq++) {
	  hPt_bq[nbq]->Fill( bq_pt[nbq] );
          hEta_bq[nbq]->Fill( bq_eta[nbq] );
          hPhi_bq[nbq]->Fill( bq_phi[nbq] );
      }

if (n_muon==4){
double muon_pt[4] = {v[muon[0]]->Pt(),v[muon[1]]->Pt(),v[muon[2]]->Pt(),v[muon[3]]->Pt()};	
        std::sort(muon_pt,muon_pt+4);
double muon_eta[4] = {v[muon[0]]->Eta(),v[muon[1]]->Eta(),v[muon[2]]->Eta(),v[muon[3]]->Eta()};	
        std::sort(muon_eta,muon_eta+4);
double muon_phi[4] = {v[muon[0]]->Phi(),v[muon[1]]->Phi(),v[muon[2]]->Phi(),v[muon[3]]->Phi()};	
        std::sort(muon_phi,muon_phi+4);

      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_muon3.SetPxPyPzE(v[muon[2]]->Px(),v[muon[2]]->Py(),v[muon[2]]->Pz(),v[muon[2]]->E());
      v_muon4.SetPxPyPzE(v[muon[3]]->Px(),v[muon[3]]->Py(),v[muon[3]]->Pz(),v[muon[3]]->E());
   
   	for (int n_m=0; n_m<4 ; n_m++) {
	  hPt_mu[n_m]->Fill( muon_pt[n_m] );
          hEta_mu[n_m]->Fill( muon_eta[n_m] );
          hPhi_mu[n_m]->Fill( muon_phi[n_m] );
      }
      v_mu_nu1.SetPxPyPzE(v[mu_nu[0]]->Px(),v[mu_nu[0]]->Py(),v[mu_nu[0]]->Pz(),v[mu_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[mu_nu[1]]->Px(),v[mu_nu[1]]->Py(),v[mu_nu[1]]->Pz(),v[mu_nu[1]]->E());
      v_mu_nu3.SetPxPyPzE(v[mu_nu[2]]->Px(),v[mu_nu[2]]->Py(),v[mu_nu[2]]->Pz(),v[mu_nu[2]]->E());
      v_mu_nu4.SetPxPyPzE(v[mu_nu[3]]->Px(),v[mu_nu[3]]->Py(),v[mu_nu[3]]->Pz(),v[mu_nu[3]]->E());


        v_w_boson1 = v_muon1+v_mu_nu1;
	v_w_boson2 = v_muon2+v_mu_nu2;
	v_w_boson3 = v_muon3+v_mu_nu3;
	v_w_boson4 = v_muon4+v_mu_nu4;

TLorentzVector v_wbosons[4]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4};
	for (int wb=0; wb<4 ; wb++) {
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
       
      }

	
	v_top1 = v_w_boson1+v_bquark1;
	v_top2 = v_w_boson2+v_bquark2;
	v_top3 = v_w_boson3+v_bquark3;
	v_top4 = v_w_boson4+v_bquark4;

}

     
if (n_muon==3){
double muon_pt[4] = {v[muon[0]]->Pt(),v[muon[1]]->Pt(),v[muon[2]]->Pt()};	
        std::sort(muon_pt,muon_pt+3);
double muon_eta[4] = {v[muon[0]]->Eta(),v[muon[1]]->Eta(),v[muon[2]]->Eta()};	
        std::sort(muon_eta,muon_eta+3);
double muon_phi[4] = {v[muon[0]]->Phi(),v[muon[1]]->Phi(),v[muon[2]]->Phi()};	
        std::sort(muon_phi,muon_phi+3);
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());
      v_muon3.SetPxPyPzE(v[muon[2]]->Px(),v[muon[2]]->Py(),v[muon[2]]->Pz(),v[muon[2]]->E());

	cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<",  "<<muon_pt[2]<<",  "<<endl;
             for (int n_m=0; n_m<3 ; n_m++) {
	  hPt_mu[n_m]->Fill( muon_pt[n_m] );
          hEta_mu[n_m]->Fill( muon_eta[n_m] );
          hPhi_mu[n_m]->Fill( muon_phi[n_m] );
          }
      v_mu_nu1.SetPxPyPzE(v[mu_nu[0]]->Px(),v[mu_nu[0]]->Py(),v[mu_nu[0]]->Pz(),v[mu_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[mu_nu[1]]->Px(),v[mu_nu[1]]->Py(),v[mu_nu[1]]->Pz(),v[mu_nu[1]]->E());
      v_mu_nu3.SetPxPyPzE(v[mu_nu[2]]->Px(),v[mu_nu[2]]->Py(),v[mu_nu[2]]->Pz(),v[mu_nu[2]]->E());
      //  w boson lorentz vectors
        v_w_boson1 = v_muon1+v_mu_nu1;
	v_w_boson2 = v_muon2+v_mu_nu2;
	v_w_boson3 = v_muon3+v_mu_nu3;


double elec_pt[4] = {v[elec[0]]->Pt()};	
        std::sort(elec_pt,elec_pt+1);
double elec_eta[4] = {v[elec[0]]->Eta()};	
        std::sort(elec_eta,elec_eta+1);
double elec_phi[4] = {v[elec[0]]->Phi()};	
        std::sort(elec_phi,elec_phi+1);
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());

      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());

      v_w_boson4 = v_elec1 + v_elec_nu1;
TLorentzVector v_wbosons[4]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4};
	for (int wb=0; wb<4 ; wb++) {
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        
    }


	v_top1 = v_w_boson1+v_bquark1;
	v_top2 = v_w_boson2+v_bquark2;
	v_top3 = v_w_boson3+v_bquark3;
	v_top4 = v_w_boson4+v_bquark4;

}

if (n_muon==2){
double muon_pt[4] = {v[muon[0]]->Pt(),v[muon[1]]->Pt()};	
        std::sort(muon_pt,muon_pt+2);
double muon_eta[4] = {v[muon[0]]->Eta(),v[muon[1]]->Eta()};	
        std::sort(muon_eta,muon_eta+2);
double muon_phi[4] = {v[muon[0]]->Phi(),v[muon[1]]->Phi()};	
        std::sort(muon_phi,muon_phi+2);

	cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<muon_pt[1]<<endl;
	for (int n_m=0; n_m<2 ; n_m++) {
	  hPt_mu[n_m]->Fill( muon_pt[n_m] );
          hEta_mu[n_m]->Fill( muon_eta[n_m] );
          hPhi_mu[n_m]->Fill( muon_phi[n_m] );
          }
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_muon2.SetPxPyPzE(v[muon[1]]->Px(),v[muon[1]]->Py(),v[muon[1]]->Pz(),v[muon[1]]->E());

      v_mu_nu1.SetPxPyPzE(v[mu_nu[0]]->Px(),v[mu_nu[0]]->Py(),v[mu_nu[0]]->Pz(),v[mu_nu[0]]->E());
      v_mu_nu2.SetPxPyPzE(v[mu_nu[1]]->Px(),v[mu_nu[1]]->Py(),v[mu_nu[1]]->Pz(),v[mu_nu[1]]->E());

  	
        v_w_boson1 = v_muon1+v_mu_nu1;
	v_w_boson2 = v_muon2+v_mu_nu2;

double elec_pt[4] = {v[elec[0]]->Pt(),v[elec[1]]->Pt()};	
        std::sort(elec_pt,elec_pt+2);
double elec_eta[4] = {v[elec[0]]->Eta(),v[elec[1]]->Eta()};	
        std::sort(elec_eta,elec_eta+2);
double elec_phi[4] = {v[elec[0]]->Phi(),v[elec[1]]->Phi()};	
        std::sort(elec_phi,elec_phi+2);

      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());

      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());

        v_w_boson3 = v_elec1 + v_elec_nu1;
        v_w_boson4 = v_elec2 + v_elec_nu2;
TLorentzVector v_wbosons[4]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4};
	for (int wb=0; wb<4 ; wb++) {
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        
    }
	
	v_top1 = v_w_boson1+v_bquark1;
	v_top2 = v_w_boson2+v_bquark2;
	v_top3 = v_w_boson3+v_bquark3;
	v_top4 = v_w_boson4+v_bquark4;

}

if (n_muon==1){
double muon_pt[4] = {v[muon[0]]->Pt()};	
        std::sort(muon_pt,muon_pt+1);
double muon_eta[4] = {v[muon[0]]->Eta()};	
        std::sort(muon_eta,muon_eta+1);
double muon_phi[4] = {v[muon[0]]->Phi()};	
        std::sort(muon_phi,muon_phi+1);

	cout<<"muon_pt :  "<<muon_pt[0]<<",  "<<endl;
	for (int n_m=0; n_m<1 ; n_m++) {
	  hPt_mu[n_m]->Fill( muon_pt[n_m] );
          hEta_mu[n_m]->Fill( muon_eta[n_m] );
          hPhi_mu[n_m]->Fill( muon_phi[n_m] );
          }
      v_muon1.SetPxPyPzE(v[muon[0]]->Px(),v[muon[0]]->Py(),v[muon[0]]->Pz(),v[muon[0]]->E());
      v_mu_nu1.SetPxPyPzE(v[mu_nu[0]]->Px(),v[mu_nu[0]]->Py(),v[mu_nu[0]]->Pz(),v[mu_nu[0]]->E());
  	
        v_w_boson1 = v_muon1+v_mu_nu1;


double elec_pt[4] = {v[elec[0]]->Pt(),v[elec[1]]->Pt(),v[elec[2]]->Pt()};	
        std::sort(elec_pt,elec_pt+3);
double elec_eta[4] = {v[elec[0]]->Eta(),v[elec[1]]->Eta(),v[elec[2]]->Eta()};	
        std::sort(elec_eta,elec_eta+3);
double elec_phi[4] = {v[elec[0]]->Phi(),v[elec[1]]->Phi(),v[elec[2]]->Phi()};	
        std::sort(elec_phi,elec_phi+3);

      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec3.SetPxPyPzE(v[elec[2]]->Px(),v[elec[2]]->Py(),v[elec[2]]->Pz(),v[elec[2]]->E());

      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      v_elec_nu3.SetPxPyPzE(v[elec_nu[2]]->Px(),v[elec_nu[2]]->Py(),v[elec_nu[2]]->Pz(),v[elec_nu[2]]->E());
       v_w_boson2 = v_elec1 + v_elec_nu1;
        v_w_boson3 = v_elec2 + v_elec_nu2;
        v_w_boson4 = v_elec3 + v_elec_nu3;
TLorentzVector v_wbosons[4]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4};
	for (int wb=0; wb<4 ; wb++) {
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        
    }
	
	v_top1 = v_w_boson1+v_bquark1;
	v_top2 = v_w_boson2+v_bquark2;
	v_top3 = v_w_boson3+v_bquark3;
	v_top4 = v_w_boson4+v_bquark4;

}

if (n_elec==4){
double elec_pt[4] = {v[elec[0]]->Pt(),v[elec[1]]->Pt(),v[elec[2]]->Pt(),v[elec[3]]->Pt()};	
        std::sort(elec_pt,elec_pt+4);
double elec_eta[4] = {v[elec[0]]->Eta(),v[elec[1]]->Eta(),v[elec[2]]->Eta(),v[elec[3]]->Eta()};	
        std::sort(elec_eta,elec_eta+4);
double elec_phi[4] = {v[elec[0]]->Phi(),v[elec[1]]->Phi(),v[elec[2]]->Phi(),v[elec[3]]->Phi()};	
        std::sort(elec_phi,elec_phi+4);
      v_elec1.SetPxPyPzE(v[elec[0]]->Px(),v[elec[0]]->Py(),v[elec[0]]->Pz(),v[elec[0]]->E());
      v_elec2.SetPxPyPzE(v[elec[1]]->Px(),v[elec[1]]->Py(),v[elec[1]]->Pz(),v[elec[1]]->E());
      v_elec3.SetPxPyPzE(v[elec[2]]->Px(),v[elec[2]]->Py(),v[elec[2]]->Pz(),v[elec[2]]->E());
      v_elec4.SetPxPyPzE(v[elec[3]]->Px(),v[elec[3]]->Py(),v[elec[3]]->Pz(),v[elec[3]]->E());

      v_elec_nu1.SetPxPyPzE(v[elec_nu[0]]->Px(),v[elec_nu[0]]->Py(),v[elec_nu[0]]->Pz(),v[elec_nu[0]]->E());
      v_elec_nu2.SetPxPyPzE(v[elec_nu[1]]->Px(),v[elec_nu[1]]->Py(),v[elec_nu[1]]->Pz(),v[elec_nu[1]]->E());
      v_elec_nu3.SetPxPyPzE(v[elec_nu[2]]->Px(),v[elec_nu[2]]->Py(),v[elec_nu[2]]->Pz(),v[elec_nu[2]]->E());
      v_elec_nu4.SetPxPyPzE(v[elec_nu[3]]->Px(),v[elec_nu[3]]->Py(),v[elec_nu[3]]->Pz(),v[elec_nu[3]]->E());

        v_w_boson1 = v_elec1 + v_elec_nu1;
        v_w_boson2 = v_elec2 + v_elec_nu2;
        v_w_boson3 = v_elec3 + v_elec_nu3;
        v_w_boson4 = v_elec4 + v_elec_nu4;
TLorentzVector v_wbosons[4]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4};
	for (int wb=0; wb<4 ; wb++) {
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        
    }
	
	v_top1 = v_w_boson1+v_bquark1;
	v_top2 = v_w_boson2+v_bquark2;
	v_top3 = v_w_boson3+v_bquark3;
	v_top4 = v_w_boson4+v_bquark4;
//cout<<"this is top quark mass:  "<<v_top4.M()<<endl;
}
//cout<<"this is outside top quark mass:  "<<v_top4.M()<<endl;

TLorentzVector v_tops[4]={v_top1, v_top2, v_top3, v_top4};
	for (int t=0; t<4 ; t++) {
        hPt_top[t]->Fill( v_tops[t].Pt() );
        hEta_top[t]->Fill( v_tops[t].Eta() );
        hPhi_top[t]->Fill( v_tops[t].Phi() );
        
    }


      int nquarks=0;
      if (q[1] != -1) nquarks = 2;
      if (q[3] != -1) nquarks = 4;


  
      for (int nq=0; nq<nquarks ; nq++) {
        hPt_q[nq]->Fill( v[q[nq]]->Pt() );
        hEta_q[nq]->Fill( v[q[nq]]->Eta() );
        hPhi_q[nq]->Fill( v[q[nq]]->Phi() );
      }

      ff>>tt;
      line++;
    
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);



    }

}
  cout << " Total number of events --> " << event << endl;


  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");

for (int i=0;i<4;i++) {
    hPt_bq[i]->Write();
    hEta_bq[i]->Write();
    hPhi_bq[i]->Write();
  }

  for (int i=0;i<4;i++) {
    hPt_q[i]->Write();
    hEta_q[i]->Write();
    hPhi_q[i]->Write();
  }

    for (int i=0;i<4;i++) {
    hPt_mu[i]->Write();
    hEta_mu[i]->Write();
    hPhi_mu[i]->Write();
  }
    for (int i=0;i<4;i++) {
    hPt_wboson[i]->Write();
    hEta_wboson[i]->Write();
    hPhi_wboson[i]->Write();
  }
    for (int i=0;i<4;i++) {
    hPt_top[i]->Write();
    hEta_top[i]->Write();
    hPhi_top[i]->Write();
  }


  rootfile->Close();

  
  exit(0);




}





  
  
