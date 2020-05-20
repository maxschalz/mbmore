#ifndef MBMORE_SRC_MULTI_ISOTOPE_ENRICH_H_
#define MBMORE_SRC_MULTI_ISOTOPE_ENRICH_H_

#include <string>

#include "cyclus.h"
#include "sim_init.h"

namespace mbmore {

class MultiIsotopeEnrich: public cyclus::Facility {
#pragma cyclus note { \
  "niche": "multi-component isotope enrichment facility", \
  "doc": "MultiIsotopeEnrich is an enrichment facility tracking minor " \
         "uranium isotopes, e.g. U234, as well. It is based on the " \
         "CYCLUS module mbmore by Meghann McGerry, Baptiste Mouginot " \
         "and Jordan Stomps." \
}
 public:
  // --- Module Members ---
  /// Constructor for the MultiIsotopeEnrich class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  MultiIsotopeEnrich(cyclus::Context* ctx);

  /// Destructor for the MultiIsotopeEnrich class
  virtual ~MultiIsotopeEnrich();

  #pragma cyclus
  /// Print information about this agent
  virtual std::string str();
  // ---

  // --- Agent Members ---
  /// Each facility is prompted to do its beginning-of-time-step
  /// stuff at the tick of the timer.
  virtual void EnterNotify();
  virtual void Tick();
  virtual void Tock();

  /// @brief The Enrichment request Materials of its given commodity.
  virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
  GetMatlRequests();

  /// @brief The Enrichment adjusts preferences for offers of
  /// natural uranium it has received to maximize U-235 content.
  /// Any offers that have zero U-235 content are not accepted.
  ///
  /// @param prefs is the preference map with all requests.
  virtual void AdjustMatlPrefs(cyclus::PrefMap<cyclus::Material>::type& prefs);

  /// @brief The Enrichment place accepted trade Materials in their Inventory.
  virtual void AcceptMatlTrades(
      const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                  cyclus::Material::Ptr> >& responses);

  /// @brief Responds to each request for this facility's commodity.  If a
  /// given request is more than this facility's inventory or SWU capacity,
  /// it will offer the minimum of its capacities.
  ///
  /// @param commod_requests is a cyclus map of requests for this facility
  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(
      cyclus::CommodMap<cyclus::Material>::type& commod_requests);

  /// @brief respond to each trade with a material enriched to the appropriate
  /// level given this facility's inventory
  ///
  /// @param trades all trades in which this trader is the supplier
  /// @param responses a container to populate with responses to each trade
  virtual void GetMatlTrades(
      const std::vector<cyclus::Trade<cyclus::Material> >& trades,
      std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                            cyclus::Material::Ptr> >& responses);
  // ---

  inline void SetMaxInventorySize(double size) {
    max_feed_inventory = size;
    inventory.capacity(size);
  }

  /// --- Convert input file flows kg/mon to SI units ---
  inline double FlowPerSec(double flow_per_mon) {
    return flow_per_mon / secpertimestep;
  }

  inline double FlowPerDt(double flow_per_sec) {
    return flow_per_sec * secpertimestep;
  }

  inline double Mg2kgPerSec(double feed_mg_per_sec) {
    return feed_mg_per_sec / (1e6);
  }
  /// ---

  /// @brief Determine if a particular material is a valid request to respond
  /// to. Valid requests must contain U235 and U238 and must have a relative
  /// U235-to-U238 ratio less than this facility's tails_assay().
  ///
  /// @return true if the above description is met by the material
  bool ValidReq(const cyclus::Material::Ptr mat);

 private:
  /// @brief Calculate the content of U235 in the feed.
  double FeedAssay(double quantity);

  /// @brief Add a material to the feed inventory
  /// @throws If the material composition is different from feed_recipe.
  void AddMat_(cyclus::Material::Ptr mat);

  /// @brief Generate a request for this facility with the quantity of
  /// required material being equal to the remaining inventory size.
  cyclus::Material::Ptr Request_();

  /// @brief Generate a material offer for a given request. The response
  /// composition will be comprised of U235, U238 and possibly of minor 
  /// uranium isotopes depending on the feed. The U235 concentration will 
  /// be as high as requested unless this enrichment is physically not 
  /// possible.
  /// @throws If the enrichment is physically not possible, e.g. producing
  /// weapons-grade from civilian spent reactor fuel.
  cyclus::Material::Ptr Offer_(cyclus::Material::Ptr req);
  
  /// @brief Perform the enrichment as requested taking into account possible
  /// feed or SWU constraints.
  cyclus::Material::Ptr Enrich_(cyclus::Material::Ptr mat, double qty);
  
  /// @brief Record the enrichment with the cyclus::Recorder
  void RecordEnrichment_(double feed_qty);

  MultiIsotopeCentrifuge centrifuge;
  MultiIsotopeCascade cascade;
  double precision = 1e-15;
  
  double sec_per_timestep = 60 * 60 * 24 * 365.25 / 12;
  
  // Set to design_tails at beginning of simulation. It gets reset if the
  // facility is used off-design.
  double tails_assay;
  
  // Set according to design assays, decimal values are rounded up.
  int n_enrich_stages;
  int n_strip_stages;
  
  #pragma cyclus var { \
    "tooltip": "feed recipe", \
    "doc": "recipe for the enrichment facility feed commodity", \
    "uilabel": "Feed Recipe", \
    "uitype": "recipe" \
  }
  std::string feed_recipe;

  #pragma cyclus var { \
    "default": 0, \
    "tooltip": "initial present uranium feed in kg", \
    "uilabel": "Initial Feed Inventory in kg", \
    "doc": "Amount of uranium feed material stored at the enrichment " \
           "facility at the beginning of the simulation in kg." \
  }
  double initial_feed;

  #pragma cyclus var { \
    "default": 1e299, \
    "tooltip": "Maximum feed material inventory in kg", \
    "uilabel": "Maximum Feed Inventory in kg", \
    "uitype": "range", \
    "range": [0.0, 1e299], \
    "doc": "Maximum total inventory of uranium feed in the enrichment " \
           "facility in kg." \
  }
  double max_feed_inventory;

  #pragma cyclus var { \
    "default": 0.000015, \
    "tooltip": "design feed flow in kg/s", \
    "uilabel": "Design Feed Flow", \
    "doc": "Target amount of feed material to be processed by the " \
    "facility in (kg/s. Either this or max_centrifuges is used to " \
    "constrain the cascade design." \
  }
  double design_feed_flow;

  #pragma cyclus var { \
    "default": -1, \
    "tooltip": "number of centrifuges available ", \
    "uilabel": "Number of Centrifuges", \
    "doc": "Number of centrifuges available to design the cascade" \
  }
  int max_centrifuges;

  #pragma cyclus var { \
    "default": 0.00711, \
    "tooltip": "design feed assay", \
    "uilabel": "Design Feed Assay", \
    "uitype": "range", \
    "range": [0.0, 1.0], \
    "doc": "Desired fraction of U235 in the feed. It should be consistent " \
           "with the feed recipe." \
  }
  double design_feed_assay;

  #pragma cyclus var { \
    "default": 0.035, \
    "tooltip": "design target product assay", \
    "uilabel": "Target Product Assay", \
    "uitype": "range", \
    "range": [0.0, 1.0], \
    "doc": "Desired fraction of U235 in the product" \
  }
  double design_product_assay;

  #pragma cyclus var { \
    "default": 0.003, \
    "tooltip": "design target tails assay", \
    "uilabel": "Target Tails Assay", \
    "uitype": "range", \
    "range": [0.0, 1.0], \
    "doc": "Desired fraction of U235 in the tails" \
  }
  double design_tails_assay;

  #pragma cyclus var { \
    "default": 320.0, \
    "tooltip": "centrifuge temperature in K", \
    "uilabel": "Centrifuge Temperature in K", \
    "doc": "Mean temperature at which the centrifuges are operated in K", \
  }
  double temperature;
  
  #pragma cyclus var { \
    "default": 485.0, \
    "tooltip": "Centrifuge velocity (m/s)", \
    "uilabel": "Centrifuge velocity (m/s)", \
    "doc": "Operational centrifuge velocity in m/s at the outer radius (a)" \
  }
  double velocity;

  #pragma cyclus var { \
    "default": 0.5, \
    "tooltip": "centrifuge height in m", \
    "uilabel": "Centrifuge Height in m", \
    "doc": "Height of the centrifuge in m" \
  }
  double height;

  #pragma cyclus var { \
    "default": 0.15, \
    "tooltip": "centrifuge diameter in m", \
    "uilabel": "Centrifuge Diameter in m", \
    "doc": "Diameter of the centrifuge in m" \
  }
  double diameter;

  #pragma cyclus var { \
    "default": 2.0, \
    "tooltip": "Centrifuge L/F* ", \
    "uilabel": "Centrifuge Countercurrent to Optimum Feed Ratio", \
    "doc": "Countercurrent to optiumum feed ratio" \
  }
  double L_over_F;

  #pragma cyclus var { \
    "default": 15.0, \
    "tooltip": "Centrifuge feed rate in mg/sec", \
    "uilabel": "Maximum feed rate for single centrifuge in mg/sec", \
    "doc" : "Maximum feed rate for a single centrifuge in mg/sec" \
  }
  double machine_feed;

  #pragma cyclus var { \
    "default": 1000.0, \
    "tooltip": "internal pressure ratio", \
    "uilabel": "Centrifuge Internal Pressure Ratio", \
    "doc": "Pressure ratio parameter for centrifuge design" \
           "that drives the r1 / r2 ratio." \
  }
  double x;

  #pragma cyclus var { \
    "default": 1.0, \
    "tooltip": "Overall centrifuge efficiency", \
    "uilabel": "Overall centrifuge efficiency", \
    "doc": "Centrifuge efficiency parameter that describes"\
           "its effective separative performance." \
  }
  double eff;

  // Input params from cycamore::Enrichment
  #pragma cyclus var { \
    "default": 1, \
    "userlevel": 10, \
    "tooltip": "Rank Material Requests by U235 Content", \
    "uilabel": "Prefer feed with higher U235 content", \
    "doc": "turn on preference ordering for input material " \
           "so that EF chooses higher U235 content first" \
  }
  bool order_prefs;

  #pragma cyclus var { \
    "default" : 1.0, \
    "tooltip" : "maximum allowed enrichment fraction", \
    "doc" : "maximum allowed weight fraction of U235 in product", \
    "uilabel" : "Maximum Allowed Enrichment", \
    "schema": '<optional>' \
              '  <element name="max_enrich">' \
              '    <data type="double">' \
              '      <param name="minInclusive">0</param>' \
              '      <param name="maxInclusive">1</param>' \
              '    </data>' \
              '  </element>' \
              '</optional>' \
  }
  double max_enrich;

  #pragma cyclus var { \
    "tooltip": "feed commodity", \
    "doc": "feed commodity that the enrichment facility accepts", \
    "uilabel": "Feed Commodity", \
    "uitype": "incommodity" \
  }
  std::string feed_commod;

  #pragma cyclus var { \
    "tooltip": "product commodity", \
    "doc": "product commodity that the enrichment facility generates", \
    "uilabel": "Product Commodity", \
    "uitype": "outcommodity" \
  }
  std::string product_commod;

  #pragma cyclus var { \
    "tooltip": "tails commodity", \
    "doc": "tails commodity supplied by enrichment facility", \
    "uilabel": "Tails Commodity", \
    "uitype": "outcommodity" \
  }
  std::string tails_commod;
  
  #pragma cyclus var{}
  cyclus::toolkit::ResBuf<cyclus::Material> inventory; // feed U
  cyclus::toolkit::ResBuf<cyclus::Material> tails; // depleted U
  
  // Used to total the intra-timestep SWU and feed U usage for meeting 
  // requests. These help to enable the time series generation.
  double intra_timestep_feed_;
  
  friend class MultiIsotopeEnrichTest;
};

} // namespace mbmore

#endif  // MBMORE_SRC_MULTI_ISOTOPE_ENRICH_H_
