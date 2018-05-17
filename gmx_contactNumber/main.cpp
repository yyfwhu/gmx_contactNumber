//
//  main.cpp
//  gmx_contactNumber
//
//  Created by Yiming Tang on 17/03/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include <iostream>
#include "contactNumber.hpp"

int main(int argc, char *argv[]) {
    
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<contactNumber>(argc, argv);
    

}
