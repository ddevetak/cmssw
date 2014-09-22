import FWCore.ParameterSet.Config as cms

PixelTrackingRecHitsValid = cms.EDAnalyzer("SiPixelTrackingRecHitsValid",
              

   src = cms.untracked.string('generalTracks'),
   runStandalone = cms.bool(True),
   outputFile = cms.untracked.string(''),
   #debugNtuple = cms.untracked.string('SiPixelTrackingRecHitsValid_Ntuple.root'),
   debugNtuple = cms.untracked.string(''),
   Fitter = cms.string('KFFittingSmoother'),
   # do we check that the simHit associated with recHit is of the expected particle type ?
   checkType = cms.bool(True),
   MTCCtrack = cms.bool(True),
   TTRHBuilder = cms.string('WithAngleAndTemplate'),


   ############################################################--------------PSet-----TH2
   
   
   TH2ResXvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                    switchon_FlippedLaddersLayer = cms.bool(False),
                                    switchon_NonFlippedLaddersLayer = cms.bool(False),
                                    switchon_total=cms.bool(False),
                                    switch_FlippedLaddersLayerTotal=cms.bool(False),
                                    switch_NonFlippedLaddersTotal=cms.bool(False),
                                    switch_LayerModule=cms.bool(False)
                                  ), 

   TH2ResYvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04),
                                    switchon_FlippedLaddersLayer = cms.bool(False),
                                    switchon_NonFlippedLaddersLayer = cms.bool(False),
                                    switchon_total=cms.bool(False),
                                    switch_FlippedLaddersLayerTotal=cms.bool(False),
                                    switch_NonFlippedLaddersTotal=cms.bool(False),
                                    switch_LayerModule=cms.bool(False)
                                  ),

   TH2ResXvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02),
                                   switchon_FlippedLaddersLayer = cms.bool(False),
                                   switchon_NonFlippedLaddersLayer = cms.bool(False),
                                   switchon_total=cms.bool(False),
                                   switch_FlippedLaddersLayerTotal=cms.bool(False),
                                   switch_NonFlippedLaddersTotal=cms.bool(False),
                                   switch_LayerModule=cms.bool(False)
                                  ),

   TH2ResYvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04),
                                   switchon_FlippedLaddersLayer = cms.bool(False),
                                   switchon_NonFlippedLaddersLayer = cms.bool(False),
                                   switchon_total=cms.bool(False),
                                   switch_FlippedLaddersLayerTotal=cms.bool(False),
                                   switch_NonFlippedLaddersTotal=cms.bool(False),
                                   switch_LayerModule=cms.bool(False)
                                  ),

   TH2PullXvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_FlippedLaddersTotal = cms.bool(True),
                                     switchon_NonFlippedLaddersTotal = cms.bool(True),
                                     switchon_LayerModule = cms.bool(True)


                                   ),

   TH2PullYvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                     switchon_total = cms.bool(True),
                                     switchon_FlippedLaddersTotal = cms.bool(True),
                                     switchon_NonFlippedLaddersTotal = cms.bool(True),
                                     switchon_LayerModule = cms.bool(True)


                                   ),

   TH2PullXvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                    switchon_total = cms.bool(True),
                                    switchon_FlippedLaddersTotal = cms.bool(True),
                                    switchon_NonFlippedLaddersTotal = cms.bool(True),
                                    switchon_LayerModule = cms.bool(True)


                                   ),

   TH2PullYvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                    switchon_total = cms.bool(True),
                                    switchon_FlippedLaddersTotal = cms.bool(True),
                                    switchon_NonFlippedLaddersTotal = cms.bool(True),
                                    switchon_LayerModule = cms.bool(True)


                                   ),


   TH2PullXvsPhiBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                   switchon_total = cms.bool(True),
                                   switchon_FlippedLaddersTotal = cms.bool(True),
                                   switchon_NonFlippedLaddersTotal = cms.bool(True),
                                   switchon_LayerModule = cms.bool(True)


                                 ),

   TH2PullYvsPhiBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                   switchon_total = cms.bool(True),
                                   switchon_FlippedLaddersTotal = cms.bool(True),
                                   switchon_NonFlippedLaddersTotal = cms.bool(True),
                                   switchon_LayerModule = cms.bool(True)


                                 ),


   TH2PullXvsEtaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(2.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                   switchon_total = cms.bool(True),
                                   switchon_FlippedLaddersTotal = cms.bool(True),
                                   switchon_NonFlippedLaddersTotal = cms.bool(True),
                                   switchon_LayerModule = cms.bool(True)


                                 ),


   TH2PullYvsEtaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(2.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                   switchon_total = cms.bool(True),
                                   switchon_FlippedLaddersTotal = cms.bool(True),
                                   switchon_NonFlippedLaddersTotal = cms.bool(True),
                                   switchon_LayerModule = cms.bool(True)


                                 ),

   TH2WPullXvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_FlippedLaddersTotal = cms.bool(False),
                                      switchon_NonFlippedLaddersTotal = cms.bool(False)

                                                  ), 

   TH2WPullYvsAlphaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(80), xmax = cms.double(100), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_FlippedLaddersTotal = cms.bool(False),
                                      switchon_NonFlippedLaddersTotal = cms.bool(False)
                                                  ), 


   TH2WPullXvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_FlippedLaddersTotal = cms.bool(False),
                                     switchon_NonFlippedLaddersTotal = cms.bool(False)
            
                                   ), 
  
   TH2WPullYvsBetaBarrel = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(10), xmax = cms.double(170), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                     switchon_FlippedLaddersTotal = cms.bool(False),
                                     switchon_NonFlippedLaddersTotal = cms.bool(False)

                                   ),

   TH2ResXvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(115.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)


                                    ),

   TH2ResYvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(115.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04),
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)


                                    ),


   TH2ResXvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67.0), xmax = cms.double(73.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)


                                   ),


   TH2ResYvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67.0), xmax = cms.double(73.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)


                                   ),

   TH2PullXvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False),
                                       switchon_DiskPlaquette = cms.bool(False)


                                     ),

   TH2PullYvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0),
                                       switchon_total = cms.bool(False),
                                       switchon_DiskPlaquette = cms.bool(False)


                                     ),

   TH2PullXvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)


                                    ),


   TH2PullYvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)


                                    ),

   TH2PullXvsPhiZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)


                                   ),

   TH2PullYvsPhiZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)


                                   ),


   TH2PullXvsEtaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(-1.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),

   TH2PullYvsEtaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(-1.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),

   TH2WPullXvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                        switchon_total = cms.bool(False)
                                       ),

   TH2WPullYvsAlphaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                        switchon_total = cms.bool(False)
                                      ),

   TH2WPullXvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False)
                                     ),

   TH2WPullYvsBetaZmPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False)

                                     ),


   TH2ResXvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(115.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)

                                    ),

   TH2ResYvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(115.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)

                                    ),

   TH2ResXvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(107), xmax = cms.double(113), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),

   TH2ResYvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(107), xmax = cms.double(113), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),

 
   TH2PullXvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False),
                                       switchon_DiskPlaquette = cms.bool(False)

                                     ),

   TH2PullYvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False),
                                       switchon_DiskPlaquette = cms.bool(False)

                                      ),


   TH2PullXvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)

                                     ),

   TH2PullYvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(False),
                                      switchon_DiskPlaquette = cms.bool(False)

                                     ),


   TH2PullXvsPhiZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),

   TH2PullYvsPhiZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),


   TH2PullXvsEtaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(1.5), xmax = cms.double(2.5), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),


   TH2PullYvsEtaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(1.5), xmax = cms.double(2.5), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(False),
                                     switchon_DiskPlaquette = cms.bool(False)

                                   ),


   TH2WPullXvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                        switchon_total = cms.bool(False)

                                      ),


   TH2WPullYvsAlphaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(100), xmax = cms.double(112.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                        switchon_total = cms.bool(False)

                                      ),


   TH2WPullXvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False)
                              
                                     ),

   TH2WPullYvsBetaZpPanel1 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(False)

                                     ),


   TH2ResXvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(65), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)


                                    ),

   TH2ResYvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(65), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)


                                    ),


   TH2ResXvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67.0), xmax = cms.double(73.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),

   TH2ResYvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67.0), xmax = cms.double(73.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2PullXvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True),
                                       switchon_DiskPlaquette = cms.bool(True)

      
                                     ),

   TH2PullYvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True),
                                       switchon_DiskPlaquette = cms.bool(True)

                                     ),


   TH2PullXvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)

                                    ),



   TH2PullYvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(72), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)

                                    ),


   TH2PullXvsPhiZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)

                                   ),


   TH2PullYvsPhiZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2PullXvsEtaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(-1.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2PullYvsEtaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-2.4), xmax = cms.double(-1.4), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),



   TH2ResXvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(65), xmax = cms.double(80.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)


                                    ),


   TH2ResYvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(65), xmax = cms.double(80.0), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)



                                    ),


   TH2ResXvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(107), xmax = cms.double(113), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.02), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2ResYvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(107), xmax = cms.double(113), Nbiny = cms.int32(100), ymin = cms.double(0.0), ymax = cms.double(0.04), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2PullXvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True),
                                       switchon_DiskPlaquette = cms.bool(True)

 
                                     ),

   TH2PullYvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80.0), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True),
                                       switchon_DiskPlaquette = cms.bool(True)



                                     ),


   TH2PullXvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)



                                    ),


   TH2PullYvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                      switchon_total = cms.bool(True),
                                      switchon_DiskPlaquette = cms.bool(True)



                                    ),


   TH2PullXvsPhiZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),



   TH2PullYvsPhiZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(-180), xmax = cms.double(180), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)



                                   ),


   TH2PullXvsEtaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(1.5), xmax = cms.double(2.5), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                   ),


   TH2PullYvsEtaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(1.5), xmax = cms.double(2.5), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True),
                                     switchon_DiskPlaquette = cms.bool(True)


                                    ),

   TH2WPullXvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True)


                                      ),

   TH2WPullYvsAlphaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                     switchon_total = cms.bool(True)


                                      ),


   TH2WPullXvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67), xmax = cms.double(73), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                     ),

   TH2WPullYvsBetaZmPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(67), xmax = cms.double(73), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                     ),


   TH2WPullXvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                      ),

   TH2WPullYvsAlphaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(68), xmax = cms.double(80), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                      ),


   TH2WPullXvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                      ),

   TH2WPullYvsBetaZpPanel2 = cms.PSet( Nbinx = cms.int32(20), xmin = cms.double(108), xmax = cms.double(112), Nbiny = cms.int32(100), ymin = cms.double(-10.0), ymax = cms.double(10.0), 
                                       switchon_total = cms.bool(True)


                                      ),





   #############################################################


   #############################################################---------------PSet-----TH1


    TH1PosxBarrel=cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1), xmax = cms.double(1), 
                            switchon_total = cms.bool(False),
                            switchon_HalfModuleTotal=cms.bool(False),
                            switchon_FullModuleTotal=cms.bool(False),
                            switchon_FlippedLaddersTotal=cms.bool(False),
                            switchon_NonFlippedLaddersTotal=cms.bool(False),
                            switchon_LayerModule=cms.bool(False)

                           ),

    TH1PosyBarrel=cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4), xmax = cms.double(4),
                            switchon_total = cms.bool(False),
                            switchon_HalfModuleTotal=cms.bool(False),
                            switchon_FullModuleTotal=cms.bool(False),
                            switchon_FlippedLaddersTotal=cms.bool(False),
                            switchon_NonFlippedLaddersTotal=cms.bool(False),
                            switchon_LayerModule=cms.bool(False)

                           ),

    TH1ErrxBarrel=cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.003),
                            switchon_total = cms.bool(False),
                            switchon_LayerModule=cms.bool(False)

                           ),

    TH1ErryBarrel=cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.01),
                            switchon_total = cms.bool(False),
                            switchon_LayerModule=cms.bool(False)

                           ),


   TH1ResxBarrel= cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.02), xmax = cms.double(0.02), 
                            switchon_Layer = cms.bool(True),
                            switchon_total=cms.bool(True),
                            switchon_LayerModule=cms.bool(True)
 
                          ),

   TH1ResyBarrel= cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.04), xmax = cms.double(0.04),
                            switchon_Layer = cms.bool(True),
                            switchon_total=cms.bool(True),
                            switchon_LayerModule=cms.bool(True)

                          ),


   TH1PullxBarrel= cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                             switchon_Layer = cms.bool(True),
                             switchon_total = cms.bool(True),
                             switchon_LayerModule = cms.bool(True)
                            
                           ),

   TH1PullyBarrel= cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0),
                             switchon_Layer = cms.bool(True),
                             switchon_total = cms.bool(True),
                             switchon_LayerModule = cms.bool(True)

                           ),


   TH1NpixBarrel = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                             switchon_total = cms.bool(False),
                             switchon_LayerModule = cms.bool(False)

                          
                           ),

   TH1NxpixBarrel = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(10.0), 
                              switchon_total = cms.bool(False),
                              switchon_LayerModule = cms.bool(False)
 
                            ),

    TH1NypixBarrel = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0),
                              switchon_total = cms.bool(False),
                              switchon_LayerModule = cms.bool(False)

                            ),
   

   TH1ChargeBarrel = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(250000.0), 
                               switchon_total = cms.bool(False),
                               switchon_LayerModule = cms.bool(False)

                             ),


   TH1PosxZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1), xmax = cms.double(1), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),

   TH1PosyZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4), xmax = cms.double(4),
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),

   TH1ErrxZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.003),
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),

  
   TH1ErryZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.01),
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),

   TH1ResxZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.02), xmax = cms.double(0.02),
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),

  
   TH1ResyZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.04), xmax = cms.double(0.04),
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),
  
   TH1PullxZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                             ),


   TH1PullyZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                             ),

  
    TH1NpixZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                             ),
 
    TH1NxpixZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(10.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                             ),


    TH1NypixZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                             ),

    TH1ChargeZmPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(100000.0),
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                                ),
   
    TH1PosxZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1), xmax = cms.double(1), 
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)
 
                               ),

    TH1PosyZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4), xmax = cms.double(4), 
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)
 
                              ),


   TH1ErrxZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.003), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),


   TH1ErryZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.01), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),


   TH1ResxZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.02), xmax = cms.double(0.02), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),


   TH1ResyZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.04), xmax = cms.double(0.04), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)

                             ),


   TH1PullxZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)
 
                              ),

   TH1PullyZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)
 
                              ),


   TH1NpixZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                               switchon_total = cms.bool(False),
                               switchon_DiskPlaquette = cms.bool(False)
 
                             ),


   TH1NxpixZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                              ),

   TH1NypixZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                                switchon_total = cms.bool(False),
                                switchon_DiskPlaquette = cms.bool(False)

                              ),


   TH1ChargeZpPanel1 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(100000.0), 
                                 switchon_total = cms.bool(False),
                                 switchon_DiskPlaquette = cms.bool(False)
 
                               ),

   #Zm-Panel2

   TH1PosxZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1), xmax = cms.double(1), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),

   TH1PosyZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4), xmax = cms.double(4), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)

                             ),


   TH1ErrxZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.003), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)
 
                             ),


   TH1ErryZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.01), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),


   TH1ResxZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.02), xmax = cms.double(0.02), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),

   TH1ResyZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.04), xmax = cms.double(0.04), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),


   TH1PullxZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)


                              ),

   TH1PullyZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)


                              ),


   TH1NpixZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                               ),


   TH1NxpixZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)

  
                              ),


   TH1NypixZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)

  
                              ),


   TH1ChargeZmPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(100000.0), 
                                 switchon_total = cms.bool(True),
                                 switchon_DiskPlaquette = cms.bool(True)


                               ),

   TH1PosxZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1), xmax = cms.double(1), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)

                              ),


   TH1PosyZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4), xmax = cms.double(4), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),

   TH1ErrxZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.003), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)

                             ),


   TH1ErryZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(0.01), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),


   TH1ResxZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.02), xmax = cms.double(0.02), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)

                             ),

   TH1ResyZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-0.04), xmax = cms.double(0.04), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),


   TH1PullxZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)

                   
                              ),


   TH1PullyZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)

                   
                              ),


   TH1NpixZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                               switchon_total = cms.bool(True),
                               switchon_DiskPlaquette = cms.bool(True)


                             ),



   TH1NxpixZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(10.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)

                   
                              ),



   TH1NypixZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(20.0), 
                                switchon_total = cms.bool(True),
                                switchon_DiskPlaquette = cms.bool(True)


                              ),


   TH1ChargeZpPanel2 = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(0.0), xmax = cms.double(100000.0), 
                                 switchon_total = cms.bool(True),
                                 switchon_DiskPlaquette = cms.bool(True)


                               ),

    
   TH1PosxBarrel_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1.0), xmax = cms.double(1.0), 
                                      switchon_total = cms.bool(True),


                                    ),

   TH1PosyBarrel_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4.0), xmax = cms.double(4.0), 
                                      switchon_total = cms.bool(True),


                                    ),

   TH1PosxZmPanel1_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1.0), xmax = cms.double(1.0), 
                                       switchon_total = cms.bool(True),


                                      ),

   TH1PosyZmPanel1_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4.0), xmax = cms.double(4.0), 
                                        switchon_total = cms.bool(True),


                                      ),


   TH1PosxZmPanel2_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1.0), xmax = cms.double(1.0), 
                                       switchon_total = cms.bool(True),


                                      ),

   TH1PosyZmPanel2_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4.0), xmax = cms.double(4.0), 
                                        switchon_total = cms.bool(True),


                                      ),

   TH1PosxZpPanel1_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1.0), xmax = cms.double(1.0),
                                       switchon_total = cms.bool(True),


                                      ),

   TH1PosyZpPanel1_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4.0), xmax = cms.double(4.0),
                                        switchon_total = cms.bool(True),


                                      ),


   TH1PosxZpPanel2_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-1.0), xmax = cms.double(1.0),
                                       switchon_total = cms.bool(True),


                                      ),

   TH1PosyZpPanel2_all_hits = cms.PSet( Nbinx = cms.int32(100), xmin = cms.double(-4.0), xmax = cms.double(4.0),
                                        switchon_total = cms.bool(True),


                                      ),


   TH1TracksPerEvent = cms.PSet( Nbinx = cms.int32(200), xmin = cms.double(0.0), xmax = cms.double(200.0),
                                        switchon_total = cms.bool(True),


                                      ),


   TH1PixRecHitsPerTrack = cms.PSet( Nbinx = cms.int32(6), xmin = cms.double(0.0), xmax = cms.double(6.0),
                                        switchon_total = cms.bool(True),


                                      ),


   #############################################################

   


                                    
   # the type of particle that the simHit associated with recHits should be
   genType = cms.int32(13),
   associatePixel = cms.bool(True),
   associateRecoTracks = cms.bool(True),
   associateStrip = cms.bool(True),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 
                         'g4SimHitsTrackerHitsPixelBarrelHighTof', 
                         'g4SimHitsTrackerHitsPixelEndcapLowTof', 
                         'g4SimHitsTrackerHitsPixelEndcapHighTof'),
   Propagator = cms.string('PropagatorWithMaterial')
                         )


