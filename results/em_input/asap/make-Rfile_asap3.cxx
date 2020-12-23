// this is the file that creates the R data object

//========================================================================
// Open the output file using the AD Model Builder template name, and
// specify 6 digits of precision
// use periods in R variable names instead of underscore

// variables used for naming fleets and indices
adstring ifleetchar;
adstring indchar; 
adstring onenum(4); 
adstring onednm(4);
adstring twodnm(4);

open_r_file(adprogram_name + ".rdat", 6, -99999);
  
  // metadata 
  open_r_info_list("info", true);
      wrt_r_item("program", "ASAP3");
  close_r_info_list();
  
    
  // basic parameter values 
  open_r_info_list("parms", false);
      wrt_r_item("styr", year1);
      wrt_r_item("endyr", (year1+nyears-1));
      wrt_r_item("nyears", nyears);
      wrt_r_item("nages", nages);
      wrt_r_item("nfleets", nfleets);
      wrt_r_item("nselblocks", nselblocks);
      wrt_r_item("navailindices", navailindices);
      wrt_r_item("nindices", nindices);
  close_r_info_list();

  // run options 
  open_r_info_list("options", false);
      wrt_r_item("isfecund", isfecund);
      wrt_r_item("frac.yr.spawn", fracyearSSB);
      wrt_r_item("do.projections", do_projections);
      wrt_r_item("ignore.guesses", ignore_guesses);
      wrt_r_item("Freport.agemin", Freport_agemin);
      wrt_r_item("Freport.agemax", Freport_agemax);
      wrt_r_item("Freport.wtopt", Freport_wtopt);
      wrt_r_item("use.likelihood.constants", use_likelihood_constants);
      wrt_r_item("Fmult.max.value", Fmult_max_value);
      wrt_r_item("N.year1.flag",NAA_year1_flag);
      wrt_r_item("do.mcmc",doMCMC);
  close_r_info_list();

  // Likelihood contributions
  open_r_info_list("like", false);
      wrt_r_item("lk.total", obj_fun);
      wrt_r_item("lk.catch.total", (lambda_catch_tot*catch_tot_likely));
      wrt_r_item("lk.discard.total", (lambda_Discard_tot*discard_tot_likely));
      wrt_r_item("lk.index.fit.total", (lambda_ind*likely_ind));
      wrt_r_item("lk.catch.age.comp", likely_catch);
      wrt_r_item("lk.discards.age.comp", likely_Discard);
      wrt_r_item("lk.index.age.comp", likely_index_age_comp);
      wrt_r_item("lk.sel.param.total", sum_sel_lambda_likely);
      wrt_r_item("lk.index.sel.param.total", sum_indexsel_lambda_likely);
      wrt_r_item("lk.q.year1", (lambda_q_year1*q_year1_likely));
      wrt_r_item("lk.q.devs", (lambda_q_devs*q_devs_likely));
      wrt_r_item("lk.Fmult.year1.total", (lambda_Fmult_year1*Fmult_year1_likely));
      wrt_r_item("lk.Fmult.devs.total", (lambda_Fmult_devs*Fmult_devs_likely));
      wrt_r_item("lk.N.year1", (lambda_N_year1_devs*N_year1_likely));
      wrt_r_item("lk.Recruit.devs", (lambda_recruit_devs*likely_SR_sigma));
      wrt_r_item("lk.SR.steepness", (lambda_steepness*steepness_likely));
      wrt_r_item("lk.SR.scaler", (lambda_SR_scaler*SR_scaler_likely));
      wrt_r_item("lk.Fmult.Max.penalty", Fmult_max_pen);
      wrt_r_item("lk.F.penalty", fpenalty);
  close_r_info_list();
 
  // fleet, block, and index specific likelihood contributions
  open_r_info_list("like.additional", false);
      wrt_r_item("nfleets",nfleets);
      wrt_r_item("nindices",nindices);
      wrt_r_item("nselparms",nselparm);
      wrt_r_item("nindexselparms",nindexselparms);
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring lk_catch_fleet = adstring("lk.catch.") + ifleetchar;
              wrt_r_item(lk_catch_fleet,(lambda_catch_tot(ifleet)*catch_tot_likely(ifleet)));
          }

          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring lk_discard_fleet = adstring("lk.discard.") + ifleetchar;
              wrt_r_item(lk_discard_fleet,(lambda_Discard_tot(ifleet)*discard_tot_likely(ifleet)));
          }

          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring lk_Fmult_year1_fleet = adstring("lk.Fmult.year1.") + ifleetchar;
              wrt_r_item(lk_Fmult_year1_fleet,(lambda_Fmult_year1(ifleet)*Fmult_year1_likely(ifleet)));
          }

          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring lk_Fmult_devs_fleet = adstring("lk.Fmult.devs.") + ifleetchar;
              wrt_r_item(lk_Fmult_devs_fleet,(lambda_Fmult_devs(ifleet)*Fmult_devs_likely(ifleet)));
          }
      }
      
      if (nindices>1)
      {
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring lk_index_fit_ind = adstring("lk.index.fit.") + indchar;
              wrt_r_item(lk_index_fit_ind,(lambda_ind(ind)*likely_ind(ind)));
          }
          
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring lk_q_year1_ind = adstring("lk.q.year1.") + indchar;
              wrt_r_item(lk_q_year1_ind,(lambda_q_year1(ind)*q_year1_likely(ind)));
          }
          
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring lk_q_devs_ind = adstring("lk.q.devs.") + indchar;
              wrt_r_item(lk_q_devs_ind,(lambda_q_devs(ind)*q_devs_likely(ind)));
          }
      }

      for (k=1;k<=nselparm;k++)
      {
          if (sel_phase(k) >=1)
          {
              if (k <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(k, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (k <=99)
              {
                  itoa(k, twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              adstring lk_sel_param = adstring("lk.sel.param.") + twodnm;
              wrt_r_item(lk_sel_param,(sel_lambda(k)*sel_likely(k)));
              
          }
      }
      
      for (k=1;k<=nindexselparms;k++)
      {
          if (indexsel_phase(k) >=1)
          {
              if (k <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(k, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (k <=99)
              {
                  itoa(k, twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              adstring lk_indexsel_param = adstring("lk.indexsel.param.") + twodnm;
              wrt_r_item(lk_indexsel_param,(indexsel_lambda(k)*indexsel_likely(k)));
              
          }
      }
      
  close_r_info_list();
  
  // initial guesses
  open_r_list("initial.guesses");
      open_r_info_list("SR.inits", false);
          wrt_r_item("is.SR.scaler.R",is_SR_scaler_R);
          wrt_r_item("SR.scaler.init",SR_scaler_ini);
          wrt_r_item("SR_steepness.init",SR_steepness_ini);
      close_r_info_list();
      wrt_r_complete_vector("NAA.year1.init",NAA_year1_ini);
      wrt_r_complete_vector("Fmult.year1.init",Fmult_year1_ini);
      wrt_r_complete_vector("q.year1.init",q_year1_ini);
      wrt_r_complete_vector("release.mort", release_mort);
      wrt_r_complete_vector("index.use.flag",use_index);
  close_r_list();
  
  // control parameters
  open_r_list("control.parms");
      open_r_info_list("phases", false);
          wrt_r_item("phase.Fmult.year1", phase_Fmult_year1);
          wrt_r_item("phase.Fmult.devs", phase_Fmult_devs);
          wrt_r_item("phase.recruit.devs", phase_recruit_devs);
          wrt_r_item("phase.N.year1.devs", phase_N_year1_devs);
          wrt_r_item("phase.q.year1", phase_q_year1);
          wrt_r_item("phase.q.devs", phase_q_devs);
          wrt_r_item("phase.SR.scaler", phase_SR_scaler);
          wrt_r_item("phase.steepness", phase_steepness);
      close_r_info_list();
      open_r_info_list("singles", false);
          wrt_r_item("lambda.N.year1.devs",lambda_N_year1_devs);
          wrt_r_item("N.year1.cv",N_year1_CV);
          wrt_r_item("lambda.recruit.devs",lambda_recruit_devs);
          wrt_r_item("lambda.steepness",lambda_steepness);
          wrt_r_item("steepness.cv",steepness_CV);
          wrt_r_item("lambda.SR.scaler",lambda_SR_scaler);
          wrt_r_item("SR.scaler.cv", SR_scaler_CV);
      close_r_info_list();
      open_r_info_list("mcmc", false);
          wrt_r_item("mcmc.nyear.opt",MCMCnyear_opt);
          wrt_r_item("mcmc.n.boot",MCMCnboot);
          wrt_r_item("mcmc.n.thin",MCMCnthin);
          wrt_r_item("mcmc.seed",MCMCseed);
          wrt_r_item("fillR.opt",fillR_opt);
          wrt_r_item("Ravg.start",Ravg_start);
          wrt_r_item("Ravg.end",Ravg_end);
      close_r_info_list();
      wrt_r_complete_vector("recruit.cv",recruit_CV);
      wrt_r_complete_vector("lambda.ind",lambda_ind);
      wrt_r_complete_vector("lambda.catch.tot",lambda_catch_tot);
      open_r_matrix("catch.tot.cv");
          wrt_r_matrix(catch_tot_CV, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nfleets);
      close_r_matrix();
      wrt_r_complete_vector("lambda.Discard.tot",lambda_Discard_tot);
      open_r_matrix("discard.tot.cv");
          wrt_r_matrix(discard_tot_CV, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nfleets);
      close_r_matrix();
      wrt_r_complete_vector("lambda.Fmult.year1",lambda_Fmult_year1);
      wrt_r_complete_vector("Fmult.year1.cv",Fmult_year1_CV);
      wrt_r_complete_vector("lambda.Fmult.devs",lambda_Fmult_devs);
      wrt_r_complete_vector("Fmult.devs.cv",Fmult_devs_CV);
      wrt_r_complete_vector("lambda.q.year1",lambda_q_year1);
      wrt_r_complete_vector("q.year1.cv",q_year1_CV);
      wrt_r_complete_vector("lambda.q.devs",lambda_q_devs);
      wrt_r_complete_vector("q.devs.cv",q_devs_CV);
      wrt_r_complete_vector("directed.fleet",directed_fleet);
      wrt_r_complete_vector("WAA.point.bio",WAApointbio);
      wrt_r_complete_vector("index.units.aggregate", index_units_aggregate);
      wrt_r_complete_vector("index.units.proportions", index_units_proportions);
      wrt_r_complete_vector("index.WAA.point", index_WAApoint);
      wrt_r_complete_vector("index.month", index_month);
      wrt_r_complete_vector("index.sel.start.age",index_start_age);
      wrt_r_complete_vector("index.sel.end.age",index_end_age);
      wrt_r_complete_vector("index.sel.choice",index_sel_choice);
      wrt_r_complete_vector("index.age.comp.flag",index_estimate_proportions);
  close_r_list();

  // selectivity input matrices for fleets and indices
  open_r_list("sel.input.mats");
      // input selectivity matrix, contains combinations of values not used, see fleet_sel_option to determine which choice was made for each block
      open_r_matrix("fleet.sel.ini");
          wrt_r_matrix(sel_ini, 2, 2);
          wrt_r_namevector(1, (nselblocks*(nages+6)));
          wrt_r_namevector(1, 4);
      close_r_matrix();
      
      open_r_matrix("index.sel.ini");
          wrt_r_matrix(index_sel_ini, 2, 2);
          wrt_r_namevector(1, (navailindices*(nages+6)));
          wrt_r_namevector(1, 4);
      close_r_matrix();
  close_r_list();
    
  // Weight at Age matrices
  open_r_list("WAA.mats");
      for (ifleet=1;ifleet<=nfleets;ifleet++)
      {
          if (nfleets < 10) itoa(ifleet, onenum, 10);
          else onenum="0";
          ifleetchar = "fleet" + onenum;
          adstring WAA_c_fleet = adstring("WAA.catch.") + ifleetchar;
          open_r_matrix(WAA_c_fleet);
              wrt_r_matrix(WAAcatchfleet(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1,nages);
          close_r_matrix();
          adstring WAA_d_fleet = adstring("WAA.discard.") + ifleetchar;
          open_r_matrix(WAA_d_fleet);
              wrt_r_matrix(WAAdiscardfleet(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1,nages);
          close_r_matrix();
      }
      open_r_matrix("WAA.catch.all");
          wrt_r_matrix(WAAcatchall, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nages);
      close_r_matrix();

      open_r_matrix("WAA.discard.all");
          wrt_r_matrix(WAAdiscardall, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nages);
      close_r_matrix();

      open_r_matrix("WAA.ssb");
          wrt_r_matrix(WAAssb, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nages);
      close_r_matrix();

      open_r_matrix("WAA.jan1");
          wrt_r_matrix(WAAjan1b, 2, 2);
          wrt_r_namevector(year1, (year1+nyears-1));
          wrt_r_namevector(1, nages);
      close_r_matrix();
      
      for (ind=1;ind<=nindices;ind++)
      {
          if (index_units_aggregate(ind)==1 || index_units_proportions(ind)==1)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring index_WAA_name = adstring("index.WAA.") + indchar;
              open_r_matrix(index_WAA_name);
                  wrt_r_matrix(index_WAA(ind), 2, 2);
                  wrt_r_namevector(year1, (year1+nyears-1));
                  wrt_r_namevector(1,nages);
              close_r_matrix();
          }
      }
      
  close_r_list();
  
  // Year by Age Matrices (not fleet specific): M, maturity, fecundity, N, Z, F, 
  open_r_matrix("M.age");
      wrt_r_matrix(M, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  open_r_matrix("maturity");
      wrt_r_matrix(mature, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  open_r_matrix("fecundity");
      wrt_r_matrix(fecundity, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  open_r_matrix("N.age");
      wrt_r_matrix(NAA, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  open_r_matrix("Z.age");
      wrt_r_matrix(Z, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  open_r_matrix("F.age");
      wrt_r_matrix(FAA_tot, 2, 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_namevector(1, nages);
  close_r_matrix();

  // Fleet by Year Matrices: Catch.tot.obs, Catch.tot.pred, Catch.tot.resid), Discard.tot.obs, Discard.tot.pred, Discard.tot.resid
  open_r_matrix("catch.obs");
      wrt_r_matrix(Catch_tot_fleet_obs, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("catch.pred");
      wrt_r_matrix(Catch_tot_fleet_pred, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("catch.std.resid");
      wrt_r_matrix(Catch_stdresid, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("discard.obs");
      wrt_r_matrix(Discard_tot_fleet_obs, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("discard.pred");
      wrt_r_matrix(Discard_tot_fleet_pred, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("discard.std.resid");
      wrt_r_matrix(Discard_stdresid, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();


  // Age Compositions: Catch and Discards observed and predicted by fleet
  open_r_list("catch.comp.mats");
      for (ifleet=1;ifleet<=nfleets;ifleet++)
      {
          if (nfleets < 10) itoa(ifleet, onenum, 10);
          else onenum="0";
          ifleetchar = "fleet" + onenum;
          adstring ccomp_ob = adstring("catch.") + ifleetchar + adstring(".ob");
          open_r_matrix(ccomp_ob);
              wrt_r_matrix(CAA_prop_obs(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(sel_start_age(ifleet), sel_end_age(ifleet));
          close_r_matrix();

          adstring ccomp_pr = adstring("catch.") + ifleetchar + adstring(".pr");
          open_r_matrix(ccomp_pr);
              wrt_r_matrix(CAA_prop_pred(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(sel_start_age(ifleet), sel_end_age(ifleet));
          close_r_matrix();

          adstring dcomp_ob = adstring("discard.") + ifleetchar + adstring(".ob");
          open_r_matrix(dcomp_ob);
              wrt_r_matrix(Discard_prop_obs(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(sel_start_age(ifleet), sel_end_age(ifleet));
          close_r_matrix();

          adstring dcomp_pr = adstring("discard.") + ifleetchar + adstring(".pr");
          open_r_matrix(dcomp_pr);
              wrt_r_matrix(Discard_prop_pred(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(sel_start_age(ifleet), sel_end_age(ifleet));
          close_r_matrix();
      }
  close_r_list();


  // fleet selectivity blocks
  open_r_matrix("fleet.sel.blocks");
      wrt_r_matrix(sel_blocks, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  // vectors of fleet selectivity options
  wrt_r_complete_vector("fleet.sel.start.age",sel_start_age);
  wrt_r_complete_vector("fleet.sel.end.age",sel_end_age);
  wrt_r_complete_vector("fleet.sel.option",sel_option);

  // selecivity matrices for each fleet 
  open_r_list("fleet.sel.mats");
      for (ifleet=1;ifleet<=nfleets;ifleet++)
      {
          if (nfleets < 10) itoa(ifleet, onenum, 10);
          else onenum="0";
          ifleetchar = "fleet" + onenum;
          adstring sel_fleet_char = adstring("sel.m.") + ifleetchar;
          open_r_matrix(sel_fleet_char);
              wrt_r_matrix(sel_by_fleet(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1, nages);
          close_r_matrix();
      }
  close_r_list();

  // Fmults by fleet
  open_r_matrix("fleet.Fmult");
      wrt_r_matrix(Fmult, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  // FAA by fleet directed and discarded
  open_r_list("fleet.FAA");
      for (ifleet=1;ifleet<=nfleets;ifleet++)
      {
          if (nfleets < 10) itoa(ifleet, onenum, 10);
          else onenum="0";
          ifleetchar = "fleet" + onenum;
          
          adstring fleet_FAA_dir = adstring("FAA.directed.") + ifleetchar;
          open_r_matrix(fleet_FAA_dir);
              wrt_r_matrix(FAA_by_fleet_dir(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1,nages);
          close_r_matrix();
          
          adstring fleet_FAA_discard = adstring("FAA.discarded.") + ifleetchar;
          open_r_matrix(fleet_FAA_discard);
              wrt_r_matrix(FAA_by_fleet_Discard(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1,nages);
          close_r_matrix();
      }
  close_r_list();
  
  // proportion release year by age matrices by fleet
  open_r_list("fleet.prop.release");
      for (ifleet=1;ifleet<=nfleets;ifleet++)
      {
          if (nfleets < 10) itoa(ifleet, onenum, 10);
          else onenum="0";
          ifleetchar = "fleet" + onenum;
          adstring fleet_prop_release = adstring("prop.release.") + ifleetchar;
          open_r_matrix(fleet_prop_release);
              wrt_r_matrix(proportion_release(ifleet), 2, 2);
              wrt_r_namevector(year1, (year1+nyears-1));
              wrt_r_namevector(1,nages);
          close_r_matrix();
      }
  close_r_list();
    
  // fleet specific annual effective sample sizes input and estimated for catch and discards
  open_r_matrix("fleet.catch.Neff.init");
      wrt_r_matrix(input_eff_samp_size_catch, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("fleet.catch.Neff.est");
      wrt_r_matrix(effective_sample_size, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("fleet.discard.Neff.init");
      wrt_r_matrix(input_eff_samp_size_discard, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  open_r_matrix("fleet.discard.Neff.est");
      wrt_r_matrix(effective_Discard_sample_size, 2, 2);
      wrt_r_namevector(1, nfleets);
      wrt_r_namevector(year1, (year1+nyears-1));
  close_r_matrix();

  // vector of q for each index if qdevs turned off, otherwise a list with vectors for each index
  if (phase_q_devs <= 0)
  {  
      wrt_r_complete_vector("q.indices",  column(q_by_index,1));
  }
  else
  {
      open_r_list("q.random.walk");
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring q_ind = adstring("q.") + indchar;
              wrt_r_complete_vector(q_ind,q_by_index(ind));
          }
      close_r_list();
  }
    
  // vectors for Freport and Biomasses (TotJan1B, SSB, ExploitableB)
  wrt_r_complete_vector("F.report",Freport);
  wrt_r_complete_vector("tot.jan1.B",TotJan1B);
  wrt_r_complete_vector("SSB",SSB);
  wrt_r_complete_vector("exploitable.B",ExploitableB);


  // F reference values 
  open_r_info_list("Fref", false);
      wrt_r_item("Fmax", Fmax_report);
      wrt_r_item("F01", F01_report);
      wrt_r_item("F30", F30SPR_report);
      wrt_r_item("F40", F40SPR_report);
      wrt_r_item("Fcurrent", Freport(nyears));
  close_r_info_list();
    
  // SR curve parameters 
  open_r_info_list("SR.parms", false);
      wrt_r_item("SR.alpha", SR_alpha);
      wrt_r_item("SR.beta", SR_beta);
      wrt_r_item("SR.SPR0", SR_spawners_per_recruit);
      wrt_r_item("SR.S0", SR_S0);
      wrt_r_item("SR.R0", SR_R0);
      wrt_r_item("SR.steepness", SR_steepness);
  close_r_info_list();

  // SR obs, pred, devs, and standardized resids
  // note year coresponds to age-1 recruitment, when plot SR curve have to offset SSB and R by one year
  open_r_df("SR.resids", year1, (year1+nyears-1), 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_df_col("year", year1, (year1+nyears-1));
      wrt_r_df_col("recruits", recruits, year1);
      wrt_r_df_col("R.no.devs", SR_pred_recruits, year1);
      wrt_r_df_col("logR.dev", log_recruit_devs, year1);
      wrt_r_df_col("SR.std.resid", SR_stdresid, year1);
  close_r_df();
    
  // annual values for S0_vec, R0_vec, steepness_vec, s_per_r_vec (last year values should match SR.parms values)
  open_r_df("SR.annual.parms", year1, (year1+nyears-1), 2);
      wrt_r_namevector(year1, (year1+nyears-1));
      wrt_r_df_col("year", year1, (year1+nyears-1));
      wrt_r_df_col("S0.vec", S0_vec, year1); 
      wrt_r_df_col("R0.vec", R0_vec, year1); 
      wrt_r_df_col("steepness.vec", steepness_vec, year1);
      wrt_r_df_col("s.per.r.vec",s_per_r_vec, year1); 
  close_r_df();



  // index stuff starts here 
  
  // selectivity by index
  open_r_matrix("index.sel");
      wrt_r_matrix(indexsel, 2, 2);
      wrt_r_namevector(1, nindices);
      wrt_r_namevector(1, nages);
  close_r_matrix();

  wrt_r_complete_vector("index.nobs",index_nobs);
    
  // index year counter (sequential numbers starting at 1 for first year)
  open_r_list("index.year.counter");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_time(ind));
      }
  close_r_list();
  
  // index years
  open_r_list("index.year");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_year(ind));
      }
  close_r_list();

  // index CV
  open_r_list("index.cv");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_cv(ind));
      }
  close_r_list();
  
  // index sigmas (derived from input CV)
  open_r_list("index.sigma");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_sigma(ind));
      }
  close_r_list();
  
  // index observations
  open_r_list("index.obs");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_obs(ind));
      }
  close_r_list();

  // predicted indices
  open_r_list("index.pred");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_pred(ind));
      }
  close_r_list();

  // index standardized residuals
  open_r_list("index.std.resid");
      for (ind=1;ind<=nindices;ind++)
      {
          if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
          {
              itoa(ind, onednm, 10);  
              twodnm = "0" + onednm;
          }
          else if (ind <=99)
          {
              itoa(ind,twodnm, 10);
          }
          else
          {
              twodnm = "00";
          }
          indchar = "ind" + twodnm;
          wrt_r_complete_vector(indchar,index_stdresid(ind));
      }
  close_r_list();
  
  // index proportions at age related output
  if (max(index_estimate_proportions)>0)  // check to see if any West Coast style indices, skip this section if all are East Coast style
  {
      // Index Age Comp   
      open_r_list("index.comp.mats");
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;

              adstring acomp_ob = indchar + adstring(".ob");
              open_r_matrix(acomp_ob);
                  wrt_r_matrix(output_index_prop_obs(ind), 2, 2);
                  wrt_r_namevector(year1, (year1+nyears-1));
                  wrt_r_namevector(1,nages);
              close_r_matrix();

              adstring acomp_pr = indchar + adstring(".pr");
              open_r_matrix(acomp_pr);
                  wrt_r_matrix(output_index_prop_pred(ind), 2, 2);
                  wrt_r_namevector(year1, (year1+nyears-1));
                  wrt_r_namevector(1, nages);
              close_r_matrix();
          }  
      close_r_list();

    // Neff for indices initial guess
    open_r_matrix("index.Neff.init");
        wrt_r_matrix(index_Neff_init, 2, 2);
        wrt_r_namevector(1, nindices);
        wrt_r_namevector(year1, (year1+nyears-1));
    close_r_matrix();

    // Neff for indices estimated
    open_r_matrix("index.Neff.est");
        wrt_r_matrix(index_Neff_est, 2, 2);
        wrt_r_namevector(1, nindices);
        wrt_r_namevector(year1, (year1+nyears-1));
    close_r_matrix();

  }  // end if-statement to test for any index age comp


  // deviations section: only reported if associated with lambda > 0
  if (lambda_N_year1_devs > 0)
  {
      // note: obs and pred include age 1 while std.resid does not - do not use age 1 when plotting
      open_r_list("deviations.N.year1");
          wrt_r_complete_vector("N.year1.obs",NAA(1));
          wrt_r_complete_vector("N.year1.pred",nyear1temp);
          wrt_r_complete_vector("N.year1.std.resid",N_year1_stdresid);
      close_r_list();
  }


  // RMSE number of observations section
  open_r_info_list("RMSE.n", false);
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_n_catch_fleet = adstring("rmse.n.catch.") + ifleetchar;
              wrt_r_item(rmse_n_catch_fleet,nyears);
          }
      }
      wrt_r_item("rmse.n.catch.tot",(nyears*nfleets));
          
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_n_discard_fleet = adstring("rmse.n.discard.") + ifleetchar;
              if (sum(Discard_tot_fleet_obs(ifleet)) > 0)
              {
                  wrt_r_item(rmse_n_discard_fleet,nyears);
              }
              else
              {
                  wrt_r_item(rmse_n_discard_fleet,0);
              }
          }
      }
      if (sum(Discard_tot_fleet_obs) > 0)
      {
          wrt_r_item("rmse.n.discard.tot",(nyears*nfleets));
      }
      else
      {
          wrt_r_item("rmse.n.discard.tot",0);
      }
          
      if (nindices>1)
      {
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring rmse_n_ind = adstring("rmse.n.") + indchar;
              wrt_r_item(rmse_n_ind,index_nobs(ind));
          }
      }
      wrt_r_item("rmse.n.ind.total",sum(index_nobs));
      
      wrt_r_item("rmse.n.N.year1",N_year1_rmse_nobs);
          
      wrt_r_item("rmse.n.Fmult.year1",Fmult_year1_rmse_nobs);
      
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_n_Fmult_devs_fleet = adstring("rmse.n.Fmult.devs.") + ifleetchar;
              wrt_r_item(rmse_n_Fmult_devs_fleet,Fmult_devs_fleet_rmse_nobs(ifleet));
          }
      }
      wrt_r_item("rmse.n.Fmult.devs.total",Fmult_devs_rmse_nobs);
      
      wrt_r_item("rmse.n.recruit.devs",SR_rmse_nobs);
      
      wrt_r_item("rmse.n.fleet.sel.params",sel_rmse_nobs);
      
      wrt_r_item("rmse.n.index.sel.params",indexsel_rmse_nobs);
      
      wrt_r_item("rmse.n.q.year1",q_year1_rmse_nobs);
      
      wrt_r_item("rmse.n.q.devs",q_devs_rmse_nobs);
      
      wrt_r_item("rmse.n.SR.steepness",steepness_rmse_nobs);
      
      wrt_r_item("rmse.n.SR.scaler",SR_scaler_rmse_nobs);
      
  close_r_info_list();
  
  // RMSE section
  open_r_info_list("RMSE", false);
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_catch_fleet = adstring("rmse.catch.") + ifleetchar;
              wrt_r_item(rmse_catch_fleet,sqrt(mean(square(Catch_stdresid(ifleet)))));
          }
      }
      wrt_r_item("rmse.catch.tot",sqrt(mean(square(Catch_stdresid))));
          
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_discard_fleet = adstring("rmse.discard.") + ifleetchar;
              if (sum(Discard_tot_fleet_obs(ifleet)) > 0)
              {
                  wrt_r_item(rmse_discard_fleet,sqrt(mean(square(Discard_stdresid(ifleet)))));
              }
              else
              {
                  wrt_r_item(rmse_discard_fleet,0);
              }
          }
      }
      if (sum(Discard_tot_fleet_obs) > 0)
      {
          wrt_r_item("rmse.discard.tot",sqrt(mean(square(Discard_stdresid))));
      }
      else
      {
          wrt_r_item("rmse.discard.tot",0);
      }
          
      if (nindices>1)
      {
          for (ind=1;ind<=nindices;ind++)
          {
              if (ind <= 9)  // note have to deal with one digit and two digit numbers separately
              {
                  itoa(ind, onednm, 10);  
                  twodnm = "0" + onednm;
              }
              else if (ind <=99)
              {
                  itoa(ind,twodnm, 10);
              }
              else
              {
                  twodnm = "00";
              }
              indchar = "ind" + twodnm;
              adstring rmse_ind = adstring("rmse.") + indchar;
              wrt_r_item(rmse_ind,sqrt(mean(square(index_stdresid(ind)))));
          }
      }
      wrt_r_item("rmse.ind.total",sqrt(mean(square(index_stdresid))));
      
      wrt_r_item("rmse.N.year1",N_year1_rmse);
          
      wrt_r_item("rmse.Fmult.year1",Fmult_year1_rmse);
      
      if (nfleets>1)
      {
          for (ifleet=1;ifleet<=nfleets;ifleet++)
          {
              if (nfleets < 10) itoa(ifleet, onenum, 10);
              else onenum="0";
              ifleetchar = "fleet" + onenum;
              adstring rmse_Fmult_devs_fleet = adstring("rmse.Fmult.devs.") + ifleetchar;
              wrt_r_item(rmse_Fmult_devs_fleet,Fmult_devs_fleet_rmse(ifleet));
          }
      }
      wrt_r_item("rmse.Fmult.devs.total",Fmult_devs_rmse);
      
      wrt_r_item("rmse.recruit.devs",SR_rmse);
      
      wrt_r_item("rmse.fleet.sel.params",sel_rmse);
      
      wrt_r_item("rmse.index.sel.params",indexsel_rmse);
      
      wrt_r_item("rmse.q.year1",q_year1_rmse);
      
      wrt_r_item("rmse.q.devs",q_devs_rmse);
      
      wrt_r_item("rmse.SR.steepness",steepness_rmse);
      
      wrt_r_item("rmse.SR.scaler",SR_scaler_rmse);
      
  close_r_info_list();
  
  open_r_list("Neff.stage2.mult");
      wrt_r_complete_vector("Neff.stage2.mult.catch", Neff_stage2_mult_catch);
      wrt_r_complete_vector("Neff.stage2.mult.discard", Neff_stage2_mult_discard);
      wrt_r_complete_vector("Neff.stage2.mult.index", Neff_stage2_mult_index);
  close_r_list();

 // close file
 close_r_file();
   
   
