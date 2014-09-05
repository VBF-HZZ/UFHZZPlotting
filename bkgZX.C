
// macro for H4l to study the Z+X background
// Author : Hengne.Li@cern.ch

{

  // Run 7, 8 or 78
  int Run = 8; // 7; 78;
 
  // integral range
  double xmin = 105.6;
  double xmax = 140.6;


  // norms 
  Double_t  nBgEvents_7TeV_4e = 1.37;
  Double_t  nBgEvents_7TeV_4mu = 0.58;
  Double_t  nBgEvents_7TeV_2e2mu = 2.29;

  Double_t  nBgEventsErr_7TeV_4e = 0.24;
  Double_t  nBgEventsErr_7TeV_4mu = 0.11;
  Double_t  nBgEventsErr_7TeV_2e2mu = 0.31;

  Double_t  nBgEvents_8TeV_4e = 6.09;
  Double_t  nBgEvents_8TeV_4mu = 3.09;
  Double_t  nBgEvents_8TeV_2e2mu = 9.25;

  Double_t  nBgEventsErr_8TeV_4e = 0.61;
  Double_t  nBgEventsErr_8TeV_4mu = 0.45;
  Double_t  nBgEventsErr_8TeV_2e2mu = 0.72;

  Double_t  nBgEvents_4e = nBgEvents_7TeV_4e+nBgEvents_8TeV_4e;
  Double_t  nBgEvents_4mu = nBgEvents_7TeV_4mu+nBgEvents_8TeV_4mu;
  Double_t  nBgEvents_2e2mu = nBgEvents_7TeV_2e2mu+nBgEvents_8TeV_2e2mu;

  Double_t  nBgEventsErr_4e = sqrt(pow(nBgEventsErr_7TeV_4e,2)+pow(nBgEventsErr_8TeV_4e,2));
  Double_t  nBgEventsErr_4mu = sqrt(pow(nBgEventsErr_7TeV_4mu,2)+pow(nBgEventsErr_8TeV_4mu,2));
  Double_t  nBgEventsErr_2e2mu = sqrt(pow(nBgEventsErr_7TeV_2e2mu,2)+pow(nBgEventsErr_8TeV_2e2mu,2));

  Double_t  nBgEvents_7TeV = nBgEvents_7TeV_4mu+nBgEvents_7TeV_4e+nBgEvents_7TeV_2e2mu;
  Double_t  nBgEvents_8TeV = nBgEvents_8TeV_4mu+nBgEvents_8TeV_4e+nBgEvents_8TeV_2e2mu;
  Double_t  nBgEvents = nBgEvents_4mu+nBgEvents_4e+nBgEvents_2e2mu;

  Double_t  nBgEventsErr_7TeV = sqrt(pow(nBgEventsErr_7TeV_4mu,2)+pow(nBgEventsErr_7TeV_4e,2)+pow(nBgEventsErr_7TeV_2e2mu,2));
  Double_t  nBgEventsErr_8TeV = sqrt(pow(nBgEventsErr_8TeV_4mu,2)+pow(nBgEventsErr_8TeV_4e,2)+pow(nBgEventsErr_8TeV_2e2mu,2));
  Double_t  nBgEventsErr = sqrt(pow(nBgEventsErr_4mu,2)+pow(nBgEventsErr_4e,2)+pow(nBgEventsErr_2e2mu,2));

  Double_t  l_4e_p[100];
  Double_t  l_4mu_p[100];
  Double_t  l_2e2mu_p[100];
  Double_t  l_4l_p[100];
  Double_t  l_4e_p_err[100];
  Double_t  l_4mu_p_err[100];
  Double_t  l_2e2mu_p_err[100];
  Double_t  l_4l_p_err[100];

  for (int i=0; i<100; i++)
  {
    l_4e_p[i] = 0.0;
    l_4mu_p[i] = 0.0;
    l_2e2mu_p[i] = 0.0;
    l_4l_p[i] = 0.0;
    l_4e_p_err[i] = 0.0;
    l_4mu_p_err[i] = 0.0;
    l_2e2mu_p_err[i] = 0.0;
    l_4l_p_err[i] = 0.0;
  }


  l_4l_p[0] = l_4e_p[0] = 1.0; //norm
  l_4l_p[1] = l_4e_p[1] = 0.117024;  // normalisation of landau1  (i.e. 2P2F)
  l_4l_p[2] = l_4e_p[2] = 195.407;   // MPV of Landau1
  l_4l_p[3] = l_4e_p[3] = 38.9472;   // sigma of Landau
  l_4l_p[4] = l_4e_p[4] = 3.68476;    // a  in pol1 = a + bx
  l_4l_p[5] = l_4e_p[5] = -0.00580439;    // b  in pol1 = a + bx
  l_4l_p[6] = l_4e_p[6] = 2.57278;     // normalisation of Landau2  (i.e. 3P1F)
  l_4l_p[7] = l_4e_p[7] = 110.862;     // MPV of Landau2
  l_4l_p[8] = l_4e_p[8] = 9.59455;     // sigma of Laudau2
  l_4l_p[9] = l_4mu_p[0] = 1.0; // norm
  l_4l_p[10] = l_4mu_p[1] = 129.0;  // MPV
  l_4l_p[11] = l_4mu_p[2] = 15.0;  // sigma
  l_4l_p[12] = l_2e2mu_p[0] = 1.0;  // norm
  l_4l_p[13] = l_2e2mu_p[1] = 0.00439895;  // normalisation of landau1  (i.e. 2P2F 2mu2e)
  l_4l_p[14] = l_2e2mu_p[2] = 195.407;     // MPV of Landau1
  l_4l_p[15] = l_2e2mu_p[3] = 38.9472;     // sigma of Landau
  l_4l_p[16] = l_2e2mu_p[4] = 3.68476;     // a  in pol1 = a + bx
  l_4l_p[17] = l_2e2mu_p[5] = -0.00580439;     // b  in pol1 = a + bx
  l_4l_p[18] = l_2e2mu_p[6] = 1.92769;     // normalisation of Landau2  (i.e. 3P1F 2mu2e)
  l_4l_p[19] = l_2e2mu_p[7] = 110.862;     // MPV of Landau2
  l_4l_p[20] = l_2e2mu_p[8] = 9.59455;     // sigma of Laudau2
  l_4l_p[21] = l_2e2mu_p[9] = 1.;               // normalisation of Landau3  (i.e. 2e2mu), set to one by convention
  l_4l_p[22] = l_2e2mu_p[10] = 129;            // MPV of Landau3
  l_4l_p[23] = l_2e2mu_p[11] = 15.0;            // sigma of Landau3

  TF1* redBgFunc_4e = new TF1("redBgFunc_4e", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6))", 10, 2000);
  for (int i=0; i<9; i++)
  {
    redBgFunc_4e->SetParameter(i,l_4e_p[i]);
    redBgFunc_4e->SetParError(i,l_4e_p_err[i]);
  } 

  TF1* redBgFunc_4mu = new TF1("redBgFunc_4mu", "landau", 10, 2000);
  for (int i=0; i<3; i++)
  {
    redBgFunc_4mu->SetParameter(i,l_4mu_p[i]);
    redBgFunc_4mu->SetParError(i,l_4mu_p_err[i]);
  }

  TF1* redBgFunc_2e2mu = new TF1("redBgFunc_2e2mu", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6) + landau(9))", 10, 2000);
  for (int i=0; i<12; i++)
  {
    redBgFunc_2e2mu->SetParameter(i,l_2e2mu_p[i]);
    redBgFunc_2e2mu->SetParError(i,l_2e2mu_p_err[i]);
  }

  //norm
  double n04e = 1./redBgFunc_4e->Integral(10, 2000);
  double n04mu = 1./redBgFunc_4mu->Integral(10, 2000);
  double n02e2mu = 1./redBgFunc_2e2mu->Integral(10, 2000);

  double n4e, n4mu, n2e2mu;
  double nerr4e, nerr4mu, nerr2e2mu;
  if (Run==7)
  {
    n4e = nBgEvents_7TeV_4e;
    n4mu = nBgEvents_7TeV_4mu;
    n2e2mu = nBgEvents_7TeV_2e2mu;
    nerr4e = nBgEventsErr_7TeV_4e;
    nerr4mu = nBgEventsErr_7TeV_4mu;
    nerr2e2mu = nBgEventsErr_7TeV_2e2mu;
  }
  else if (Run==8)
  {
    n4e = nBgEvents_8TeV_4e;
    n4mu = nBgEvents_8TeV_4mu;
    n2e2mu = nBgEvents_8TeV_2e2mu;
    nerr4e = nBgEventsErr_8TeV_4e;
    nerr4mu = nBgEventsErr_8TeV_4mu;
    nerr2e2mu = nBgEventsErr_8TeV_2e2mu;
  }
  else if (Run==78)
  {
    n4e = nBgEvents_4e;
    n4mu = nBgEvents_4mu;
    n2e2mu = nBgEvents_2e2mu;
    nerr4e = nBgEventsErr_4e;
    nerr4mu = nBgEventsErr_4mu;
    nerr2e2mu = nBgEventsErr_2e2mu;
  }

  l_4l_p[0] = l_4e_p[0] = n4e*n04e;
  l_4l_p[9] = l_4mu_p[0] = n4mu*n04mu;
  l_4l_p[12] = l_2e2mu_p[0] = n2e2mu*n02e2mu;

  l_4l_p_err[0] = l_4e_p_err[0] = nerr4e*n04e;
  l_4l_p_err[9] = l_4mu_p_err[0] = nerr4mu*n04mu;
  l_4l_p_err[12] = l_2e2mu_p_err[0] = nerr2e2mu*n02e2mu;

  redBgFunc_4e->SetParameter(0,l_4e_p[0]);
  redBgFunc_4mu->SetParameter(0,l_4mu_p[0]);
  redBgFunc_2e2mu->SetParameter(0,l_2e2mu_p[0]);
  redBgFunc_4e->SetParError(0, l_4e_p_err[0]);
  redBgFunc_4mu->SetParError(0, l_4mu_p_err[0]);
  redBgFunc_2e2mu->SetParError(0, l_2e2mu_p_err[0]);

  // 4l
  TF1* redBgFunc_4l = new TF1("redBgFunc_4l", "[0]*(landau(1) * (1 + exp( pol1(4))) + landau(6))+landau(9)+[12]*(landau(13) * (1 + exp( pol1(16))) + landau(18) + landau(21))", 10, 2000);
  for (int i=0; i<24; i++)
  {
    redBgFunc_4l->SetParameter(i,l_4l_p[i]);
    redBgFunc_4l->SetParError(i,l_4l_p_err[i]);
  }

  redBgFunc_4l->SetRange(xmin,xmax);
  redBgFunc_4l->Draw();
  
  double Nint_4e, Nint_4mu, Nint_2e2mu, Nint_4l;
  double NintErr_4e, NintErr_4mu, NintErr_2e2mu, NintErr_4l;


  // making cov matrix
  Double_t l_4e_cov[9*9];
  Double_t l_4mu_cov[3*3];
  Double_t l_2e2mu_cov[12*12];
  Double_t l_4l_cov[24*24];

  for (int i=0; i<9; i++) for (int j=0; j<9; j++) l_4e_cov[i*9+j] = l_4e_p_err[i]*l_4e_p_err[j];
  for (int i=0; i<3; i++) for (int j=0; j<3; j++) l_4mu_cov[i*3+j] = l_4mu_p_err[i]*l_4mu_p_err[j];
  for (int i=0; i<12; i++) for (int j=0; j<12; j++) l_2e2mu_cov[i*12+j] = l_2e2mu_p_err[i]*l_2e2mu_p_err[j];
  for (int i=0; i<24; i++) for (int j=0; j<24; j++) l_4l_cov[i*24+j] = l_4l_p_err[i]*l_4l_p_err[j];

  // get integrals and integral errors
  Nint_4e = redBgFunc_4e->Integral(xmin, xmax);
  Nint_4mu = redBgFunc_4mu->Integral(xmin, xmax);
  Nint_2e2mu = redBgFunc_2e2mu->Integral(xmin, xmax);
  Nint_4l = redBgFunc_4l->Integral(xmin, xmax);
  NintErr_4e = redBgFunc_4e->IntegralError(xmin, xmax, l_4e_p, l_4e_cov);
  NintErr_4mu = redBgFunc_4mu->IntegralError(xmin, xmax, l_4mu_p, l_4mu_cov);
  NintErr_2e2mu = redBgFunc_2e2mu->IntegralError(xmin, xmax, l_2e2mu_p, l_2e2mu_cov);
  NintErr_4l = redBgFunc_4l->IntegralError(xmin, xmax, l_4l_p, l_4l_cov);

  // print
  std::cout << " Z+X background for " << Run << " TeV in range " << xmin << " < m4l < " << xmax << ":" << std::endl;
  std::cout << " N_4e    : " << Nint_4e    << " +- " << NintErr_4e    << std::endl;
  std::cout << " N_4mu   : " << Nint_4mu   << " +- " << NintErr_4mu   << std::endl;
  std::cout << " N_2e2mu : " << Nint_2e2mu << " +- " << NintErr_2e2mu << std::endl;
  std::cout << " N_4l    : " << Nint_4l    << " +- " << NintErr_4l    << std::endl;
}
