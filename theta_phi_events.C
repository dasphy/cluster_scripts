#define theta_phi_events_cxx
#include "theta_phi_events.h"
#include <TH2.h>
#include <TStyle.h>

void theta_phi_events::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void theta_phi_events::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}

Bool_t theta_phi_events::Process(Long64_t entry)
{
   fReader.SetLocalEntry(entry);
   //std::cout << "Size is: " << CaloClusters_position_x.GetSize() << std::endl;
   double particle_theta = atan2(sqrt(genParticles_momentum_y[0]*genParticles_momentum_y[0]+genParticles_momentum_x[0]*genParticles_momentum_x[0]), genParticles_momentum_z[0]);
   double particle_phi = atan2(genParticles_momentum_y[0], genParticles_momentum_x[0]);
   P_theta->Fill(particle_theta);
   P_phi->Fill(particle_phi);

   for (int i=0; i<6; i++) {
     // 6 edges of 5 sub-ranges
     theta_range.push_back(M_PI/180.*40. + i * M_PI/180.*(140.-40.)/5.);
     phi_range.push_back(-M_PI + i * 2.*M_PI/5.);
     if ( i<5 ) {
       // middle of sub-ranges, for plotting TGraph
       x_axis_theta.push_back(M_PI/180.*40. + (i+.5) * M_PI/180.*(140.-40.)/5.);
       x_axis_phi.push_back(-M_PI + (i+.5) * 2.*M_PI/5.);
     }
   }

   for (int i=0; i<5; i++) {
     if ( fixed_phi ? (particle_theta > theta_range[i] && particle_theta < theta_range[i+1]) : (particle_phi > phi_range[i] && particle_phi < phi_range[i+1]) ) {
       for (unsigned int j = 0; j < CaloClusters_position_x.GetSize(); j ++) {
         double cluster_theta = atan2(sqrt(CaloClusters_position_y[j]*CaloClusters_position_y[j]+CaloClusters_position_x[j]*CaloClusters_position_x[j]), CaloClusters_position_z[j]);
         double cluster_phi = atan2(CaloClusters_position_y[j], CaloClusters_position_x[j]);

         if (i==0) {
           C_theta_1->Fill(cluster_theta - particle_theta);
           C_phi_1->Fill(cluster_phi - particle_phi);
         }
         if (i==1) {
           C_theta_2->Fill(cluster_theta - particle_theta);
           C_phi_2->Fill(cluster_phi - particle_phi);
         }
         if (i==2) {
           C_theta_3->Fill(cluster_theta - particle_theta);
           C_phi_3->Fill(cluster_phi - particle_phi);
         }
         if (i==3) {
           C_theta_4->Fill(cluster_theta - particle_theta);
           C_phi_4->Fill(cluster_phi - particle_phi);
         }
         if (i==4) {
           C_theta_5->Fill(cluster_theta - particle_theta);
           C_phi_5->Fill(cluster_phi - particle_phi);
         }
       }
     }
   }

   return kTRUE;
}

void theta_phi_events::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void theta_phi_events::Terminate()
{
   //double C_theta_GetStdDev = C_theta->GetStdDev();
   //std::cout << "C_theta_GetStdDev == " << C_theta_GetStdDev << std::endl;

   if (fixed_phi == false) {
     label = "fix_theta";
   } else {
     label = "fix_phi";
   }

   TCanvas* c = new TCanvas("canvas", "canvas", 800, 600);

   P_theta->Draw();
   c->Print("./plot/particle_theta_"+ label +".png");
   P_phi->Draw();
   c->Print("./plot/particle_phi_"+ label +".png");

   C_theta_1->Fit("gaus");
   C_theta_1->Draw();
   c->Print("./plot/response_theta_1_"+ label +".png");
   C_theta_2->Fit("gaus");
   C_theta_2->Draw();
   c->Print("./plot/response_theta_2_"+ label +".png");
   C_theta_3->Fit("gaus");
   C_theta_3->Draw();
   c->Print("./plot/response_theta_3_"+ label +".png");
   C_theta_4->Fit("gaus");
   C_theta_4->Draw();
   c->Print("./plot/response_theta_4_"+ label +".png");
   C_theta_5->Fit("gaus");
   C_theta_5->Draw();
   c->Print("./plot/response_theta_5_"+ label +".png");

   C_phi_1->Fit("gaus");
   C_phi_1->Draw();
   c->Print("./plot/response_phi_1_"+ label +".png");
   C_phi_2->Fit("gaus");
   C_phi_2->Draw();
   c->Print("./plot/response_phi_2_"+ label +".png");
   C_phi_3->Fit("gaus");
   C_phi_3->Draw();
   c->Print("./plot/response_phi_3_"+ label +".png");
   C_phi_4->Fit("gaus");
   C_phi_4->Draw();
   c->Print("./plot/response_phi_4_"+ label +".png");
   C_phi_5->Fit("gaus");
   C_phi_5->Draw();
   c->Print("./plot/response_phi_5_"+ label +".png");

   TF1 *tf_C_theta_1 = (TF1*)C_theta_1->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_theta_1 = tf_C_theta_1->GetParameter(2);
   resol_theta.push_back(sigma_C_theta_1);
   TF1 *tf_C_theta_2 = (TF1*)C_theta_2->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_theta_2 = tf_C_theta_2->GetParameter(2);
   resol_theta.push_back(sigma_C_theta_2);
   TF1 *tf_C_theta_3 = (TF1*)C_theta_3->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_theta_3 = tf_C_theta_3->GetParameter(2);
   resol_theta.push_back(sigma_C_theta_3);
   TF1 *tf_C_theta_4 = (TF1*)C_theta_4->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_theta_4 = tf_C_theta_4->GetParameter(2);
   resol_theta.push_back(sigma_C_theta_4);
   TF1 *tf_C_theta_5 = (TF1*)C_theta_5->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_theta_5 = tf_C_theta_5->GetParameter(2);
   resol_theta.push_back(sigma_C_theta_5);

   TF1 *tf_C_phi_1 = (TF1*)C_phi_1->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_phi_1 = tf_C_phi_1->GetParameter(2);
   resol_phi.push_back(sigma_C_phi_1);
   TF1 *tf_C_phi_2 = (TF1*)C_phi_2->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_phi_2 = tf_C_phi_2->GetParameter(2);
   resol_phi.push_back(sigma_C_phi_2);
   TF1 *tf_C_phi_3 = (TF1*)C_phi_3->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_phi_3 = tf_C_phi_3->GetParameter(2);
   resol_phi.push_back(sigma_C_phi_3);
   TF1 *tf_C_phi_4 = (TF1*)C_phi_4->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_phi_4 = tf_C_phi_4->GetParameter(2);
   resol_phi.push_back(sigma_C_phi_4);
   TF1 *tf_C_phi_5 = (TF1*)C_phi_5->GetListOfFunctions()->FindObject("gaus");
   double sigma_C_phi_5 = tf_C_phi_5->GetParameter(2);
   resol_phi.push_back(sigma_C_phi_5);

   std::vector<double> error_x_theta =  {M_PI/180.*(140.-40.)/10., M_PI/180.*(140.-40.)/10., M_PI/180.*(140.-40.)/10., M_PI/180.*(140.-40.)/10., M_PI/180.*(140.-40.)/10.};
   std::vector<double> error_x_phi =  {M_PI/5., M_PI/5., M_PI/5., M_PI/5., M_PI/5.};
   std::vector<double> error_y =  {0., 0., 0., 0., 0.};

   TCanvas* c2 = new TCanvas("canvas2", "canvas2", 800, 600);

   TGraphErrors gr_resol_theta(5, &x_axis_theta[0], &resol_theta[0], &error_x_theta[0], &error_y[0]);
   TGraphErrors gr_resol_phi(5, &x_axis_phi[0], &resol_phi[0], &error_x_phi[0], &error_y[0]);
   gr_resol_theta.SetMarkerStyle(34);
   gr_resol_theta.SetMarkerColor(kRed);
   gr_resol_phi.SetMarkerStyle(34);
   gr_resol_phi.SetMarkerColor(kRed);

   TH2D* axor1 = new TH2D("axor1", "", 1000, M_PI/180.*40., M_PI/180.*140., 1000, 0., 2.*resol_theta[2]);
   axor1->SetDirectory(0);
   axor1->GetXaxis()->SetTitle("#theta");
   axor1->GetYaxis()->SetTitle("Resolution (#sigma)");
   //axor1->GetYaxis()->SetTitleOffset(0.5);
   axor1->SetStats(0);
   axor1->Draw();
   gr_resol_theta.DrawClone("SAME P");
   c2->Print("./plot/graph_resol_theta_"+ label +".png");

   TH2D* axor2 = new TH2D("axor2", "", 1000, -M_PI, M_PI, 1000, 0., 2.*resol_phi[2]);
   axor2->SetDirectory(0);
   axor2->GetXaxis()->SetTitle("#phi");
   axor2->GetYaxis()->SetTitle("Resolution (#sigma)");
   axor2->SetStats(0);
   axor2->Draw();
   gr_resol_phi.DrawClone("SAME P");
   c2->Print("./plot/graph_resol_phi_"+ label +".png");
}
