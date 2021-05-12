//
// Created by Andrea Bonvini on 11/05/21.
//

#ifndef POINTPROCESS_SERIALIZE_H
#define POINTPROCESS_SERIALIZE_H


void ppRes2csv(PointProcessResult& ppRes, const std::string& outputResultsName, const std::string& outputTausName){
    // The main Regression Results will be saved at outputResultsName
    // The taus will be saved in outputTausName
    std::ofstream csv(outputResultsName);
    std::ofstream taus(outputTausName);
    // Define header...
    csv << ((ppRes.distribution == PointProcessDistributions::InverseGaussian) ? "Time,Mu,Sigma,Kappa,Lambda,eventHappened,nIter,Theta0" : "Time,Mu,Sigma,Lambda,eventHappened,nIter,Theta0");
    for(long i = 0 ; i < ppRes.AR_ORDER; i++) {
        csv << "," << "Theta" << std::to_string(i + 1);
    }
    csv << "\n";

    if (ppRes.distribution == PointProcessDistributions::InverseGaussian) {
        // In case we are serializing an InverseGaussian Regression Result we have to save also the Kappa parameter for
        // each time step.
        for (auto &res: ppRes.results) {
            auto tmp = dynamic_cast<IGRegressionResult *>(res.get());
            csv << ppRes.t0 + tmp->time << "," << tmp->mu << "," << tmp->sigma << "," << tmp->kappa << "," << tmp->lambda << ","
                << tmp->eventHappened << "," << tmp->nIter << "," << tmp->theta0;
            for (long i = 0; i < tmp->thetaP.size(); i++) {
                csv << "," << tmp->thetaP(i);
            }
            csv << "\n";
        }
    }
    else{
        for (auto& res: ppRes.results){
            csv << ppRes.t0 + res->time << "," << res->mu << "," << res->sigma << "," << res->lambda << "," << res->eventHappened << "," << res->nIter << "," << res->theta0 ;
            for (long i = 0 ; i < res->thetaP.size(); i++){
                csv << "," << res->thetaP(i);
            }
            csv << "\n";
        }
    }

    taus << "Taus\n";
    for (auto& tau: ppRes.taus){
        taus<< tau << "\n";
    }
    csv.close();
    taus.close();

}

#endif //POINTPROCESS_SERIALIZE_H
