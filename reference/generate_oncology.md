# Generate Dataset that simulates an oncology trial

Generate Dataset that simulates an oncology trial

Create an empty assumptions data.frame for generate_oncology

Calculate true summary statistics for scenarios with delayed treatment
effect

## Usage

``` r
generate_oncology(
  condition,
  fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL,
    allow_random_cens = TRUE)
)

oncology_scenario()

oncology_scenario_set_truevalues(Design, fixed_objects = NULL)
```

## Arguments

- condition:

  condition row of Design dataset

- fixed_objects:

  fixed objects not used for now

- Design:

  Design data.frame for x

## Value

For generate_oncology: A data set with n rows and the columns .....

For oncology_scenario: a design tibble with default values invisibly

For

## Details

Condition has to contain the following columns:

- k \# number of post baseline visits

- recr_interval \# recruitment interval in years (e.g. 2 means that
  recruitment takes 2 years)

- max_duration \# maximal duration of the trial in years (e.g. 7 means
  that the last patient recruited is followed for 5 years)

- aimed_for_n_per_required_event \# This governs the recruitment rate.

- alpha \# type I error rate for the logrank test at the end of the
  trial (i.e. at max_duration)

- power \# power for the logrank test at the end of the trial (i.e. at
  max_duration)

- p_trt \# proportion of patients randomized to the treatment arm

- w \# threshold for covariate W to determine the effect of W on the
  hazard

- mu_W \# list with list entries containing mean values of W at each
  visit for control and treatment group in k dimensional vectors

- mu_L \# list with list entries containing mean values of L at each
  visit for control and treatment group in k dimensional vectors

- Sigma_W_L \# list with matrix of dimension k x k with the covariance
  matrix of W and L at all visits (assuming the same covariance matrix
  for control and treatment group)

- beta_prog \# list with vectors of dimension 7 with the coefficients
  for the progression hazard at each visit for the variables Int, X, W,
  W\>0, L, trt, switched

- beta_switch \# list with vectors of dimension 7 with the coefficients
  for the switching probability at each visit for the variables Int, X,
  W, W\>0, L, trt, switched

- beta_cens \# list with vectors of dimension 7 with the coefficients
  for the censoring hazard at each visit for the variables Int, X, W,
  W\>0, L, trt, switched

- beta_death \# list with vectors of dimension 8 with the coefficients
  for the death hazard at each visit for the variables Int, X, W, W\>0,
  L, trt, switched, logHR_assumed

oncology_scenario generates a default design `data.frame` for use with
generate_oncology.

## Functions

- `generate_oncology()`: simulates a data set with n rows.

- `oncology_scenario()`: generate default design tibble

- `oncology_scenario_set_truevalues()`: calculate true summary
  statistics for ...

## Examples

``` r
Design <- oncology_scenario()
generate_oncology(Design[1, ])
#>      id trt          X_0          X_1          X_2         X_3 X_4 X_5 X_6 X_7
#> 73    1   0 -1.327331255           NA           NA          NA  NA  NA  NA  NA
#> 89    2   0 -0.265020949           NA           NA          NA  NA  NA  NA  NA
#> 71    3   0 -1.129568395           NA           NA          NA  NA  NA  NA  NA
#> 120   4   1  0.792347143           NA           NA          NA  NA  NA  NA  NA
#> 77    5   0 -2.128069870           NA           NA          NA  NA  NA  NA  NA
#> 54    6   0  0.485275870           NA           NA          NA  NA  NA  NA  NA
#> 121   7   1 -1.177596300           NA           NA          NA  NA  NA  NA  NA
#> 141   8   1 -0.762932289           NA           NA          NA  NA  NA  NA  NA
#> 70    9   0  0.048068709           NA           NA          NA  NA  NA  NA  NA
#> 34   10   0 -0.036213207           NA           NA          NA  NA  NA  NA  NA
#> 50   11   0 -1.624615430           NA           NA          NA  NA  NA  NA  NA
#> 111  12   1 -1.307237931           NA           NA          NA  NA  NA  NA  NA
#> 94   13   0 -0.734696342           NA           NA          NA  NA  NA  NA  NA
#> 107  14   1 -0.725479323           NA           NA          NA  NA  NA  NA  NA
#> 19   15   0 -0.775199253           NA           NA          NA  NA  NA  NA  NA
#> 115  16   1 -1.230762282           NA           NA          NA  NA  NA  NA  NA
#> 23   17   0 -0.371732264 -0.371732264           NA          NA  NA  NA  NA  NA
#> 47   18   0 -1.543695084           NA           NA          NA  NA  NA  NA  NA
#> 67   19   0 -0.417580500           NA           NA          NA  NA  NA  NA  NA
#> 81   20   0 -0.115655250           NA           NA          NA  NA  NA  NA  NA
#> 142  21   1 -0.464005889           NA           NA          NA  NA  NA  NA  NA
#> 37   22   0 -0.347298344           NA           NA          NA  NA  NA  NA  NA
#> 8    23   0 -1.518914931           NA           NA          NA  NA  NA  NA  NA
#> 79   24   0 -2.299280635 -2.299280635           NA          NA  NA  NA  NA  NA
#> 95   25   0  0.814690644  0.814690644           NA          NA  NA  NA  NA  NA
#> 26   26   0 -0.595717879           NA           NA          NA  NA  NA  NA  NA
#> 155  27   1 -2.149099081 -2.149099081           NA          NA  NA  NA  NA  NA
#> 48   28   0 -0.811621646           NA           NA          NA  NA  NA  NA  NA
#> 158  29   1 -0.961820003           NA           NA          NA  NA  NA  NA  NA
#> 200  30   1 -1.096651989           NA           NA          NA  NA  NA  NA  NA
#> 134  31   1 -0.845100201           NA           NA          NA  NA  NA  NA  NA
#> 156  32   1 -0.169367187           NA           NA          NA  NA  NA  NA  NA
#> 119  33   1 -0.693481350           NA           NA          NA  NA  NA  NA  NA
#> 127  34   1 -1.775397920 -1.775397920           NA          NA  NA  NA  NA  NA
#> 103  35   1  0.415851412           NA           NA          NA  NA  NA  NA  NA
#> 193  36   1  0.210213656  0.210213656           NA          NA  NA  NA  NA  NA
#> 118  37   1  0.132161461           NA           NA          NA  NA  NA  NA  NA
#> 209  38   1 -1.534128183           NA           NA          NA  NA  NA  NA  NA
#> 98   39   0 -0.520123482           NA           NA          NA  NA  NA  NA  NA
#> 61   40   0 -1.793427835           NA           NA          NA  NA  NA  NA  NA
#> 4    41   0  1.920388895           NA           NA          NA  NA  NA  NA  NA
#> 157  42   1 -0.829284774           NA           NA          NA  NA  NA  NA  NA
#> 7    43   0 -0.803455568           NA           NA          NA  NA  NA  NA  NA
#> 66   44   0 -2.696984957           NA           NA          NA  NA  NA  NA  NA
#> 60   45   0 -0.657858374 -0.657858374           NA          NA  NA  NA  NA  NA
#> 43   46   0  0.226715318  0.226715318           NA          NA  NA  NA  NA  NA
#> 210  47   1  1.936517565           NA           NA          NA  NA  NA  NA  NA
#> 28   48   0  0.677246175  0.677246175           NA          NA  NA  NA  NA  NA
#> 17   49   0 -0.711604050 -0.711604050           NA          NA  NA  NA  NA  NA
#> 162  50   1 -0.303133823           NA           NA          NA  NA  NA  NA  NA
#> 205  51   1  0.743963173           NA           NA          NA  NA  NA  NA  NA
#> 80   52   0 -0.631220028           NA           NA          NA  NA  NA  NA  NA
#> 188  53   1  0.860846390  0.860846390           NA          NA  NA  NA  NA  NA
#> 25   54   0  0.039356836           NA           NA          NA  NA  NA  NA  NA
#> 84   55   0 -1.291981751 -1.291981751           NA          NA  NA  NA  NA  NA
#> 108  56   1 -0.523347405 -0.523347405 -0.523347405          NA  NA  NA  NA  NA
#> 9    57   0 -1.555776031           NA           NA          NA  NA  NA  NA  NA
#> 144  58   1 -0.213436234           NA           NA          NA  NA  NA  NA  NA
#> 195  59   1  0.830192409  0.830192409           NA          NA  NA  NA  NA  NA
#> 168  60   1  0.126316649  0.126316649           NA          NA  NA  NA  NA  NA
#> 149  61   1 -0.671788473 -0.671788473 -0.671788473          NA  NA  NA  NA  NA
#> 41   62   0 -0.775025328 -0.775025328           NA          NA  NA  NA  NA  NA
#> 113  63   1  0.164531916  0.164531916           NA          NA  NA  NA  NA  NA
#> 167  64   1 -0.325262792 -0.325262792 -0.325262792          NA  NA  NA  NA  NA
#> 44   65   0 -0.711380049 -0.711380049           NA          NA  NA  NA  NA  NA
#> 166  66   1 -1.425972572 -1.425972572           NA          NA  NA  NA  NA  NA
#> 55   67   0 -0.132628455 -0.132628455           NA          NA  NA  NA  NA  NA
#> 92   68   0  0.532975205  0.532975205           NA          NA  NA  NA  NA  NA
#> 197  69   1  0.698499194  0.698499194  0.698499194          NA  NA  NA  NA  NA
#> 143  70   1 -1.950884497 -1.950884497 -1.950884497          NA  NA  NA  NA  NA
#> 12   71   0 -1.343792501 -1.343792501           NA          NA  NA  NA  NA  NA
#> 164  72   1  0.301795118  0.301795118  0.301795118          NA  NA  NA  NA  NA
#> 22   73   0 -0.299127703 -0.299127703           NA          NA  NA  NA  NA  NA
#> 42   74   0 -0.153702393 -0.153702393 -0.153702393          NA  NA  NA  NA  NA
#> 46   75   0  0.084716934  0.084716934           NA          NA  NA  NA  NA  NA
#> 68   76   0  0.175446188  0.175446188           NA          NA  NA  NA  NA  NA
#> 18   77   0  1.005909127  1.005909127           NA          NA  NA  NA  NA  NA
#> 179  78   1  0.104606714  0.104606714  0.104606714          NA  NA  NA  NA  NA
#> 75   79   0  1.086098823  1.086098823  1.086098823          NA  NA  NA  NA  NA
#> 96   80   0  0.629507845  0.629507845  0.629507845  0.62950785  NA  NA  NA  NA
#> 49   81   0 -0.223436252 -0.223436252           NA          NA  NA  NA  NA  NA
#> 33   82   0 -0.232978661 -0.232978661 -0.232978661          NA  NA  NA  NA  NA
#> 129  83   1 -1.208932144 -1.208932144           NA          NA  NA  NA  NA  NA
#> 122  84   1  0.194159505  0.194159505           NA          NA  NA  NA  NA  NA
#> 102  85   1 -0.704295442 -0.704295442 -0.704295442 -0.70429544  NA  NA  NA  NA
#> 51   86   0 -1.298458855 -1.298458855 -1.298458855 -1.29845886  NA  NA  NA  NA
#> 88   87   0  0.823198729  0.823198729           NA          NA  NA  NA  NA  NA
#> 161  88   1  0.521056693  0.521056693           NA          NA  NA  NA  NA  NA
#> 160  89   1 -0.784471147 -0.784471147 -0.784471147          NA  NA  NA  NA  NA
#> 105  90   1  0.898586893  0.898586893           NA          NA  NA  NA  NA  NA
#> 86   91   0 -1.189930229 -1.189930229           NA          NA  NA  NA  NA  NA
#> 10   92   0 -0.534201600 -0.534201600 -0.534201600 -0.53420160  NA  NA  NA  NA
#> 15   93   0 -1.790069456 -1.790069456 -1.790069456          NA  NA  NA  NA  NA
#> 93   94   0 -0.202056466 -0.202056466 -0.202056466 -0.20205647  NA  NA  NA  NA
#> 64   95   0  1.308854820  1.308854820  1.308854820          NA  NA  NA  NA  NA
#> 110  96   1 -1.744183492 -1.744183492 -1.744183492 -1.74418349  NA  NA  NA  NA
#> 62   97   0 -1.033135468 -1.033135468           NA          NA  NA  NA  NA  NA
#> 172  98   1  0.290301415  0.290301415  0.290301415  0.29030141  NA  NA  NA  NA
#> 174  99   1 -0.834203453 -0.834203453 -0.834203453          NA  NA  NA  NA  NA
#> 192 100   1  0.572901813  0.572901813  0.572901813          NA  NA  NA  NA  NA
#> 106 101   1  1.332554213  1.332554213  1.332554213          NA  NA  NA  NA  NA
#> 82  102   0  0.690667743  0.690667743  0.690667743  0.69066774  NA  NA  NA  NA
#> 190 103   1 -1.745148998 -1.745148998 -1.745148998 -1.74514900  NA  NA  NA  NA
#> 39  104   0 -1.151406159 -1.151406159 -1.151406159 -1.15140616  NA  NA  NA  NA
#> 45  105   0 -0.473650038 -0.473650038 -0.473650038 -0.47365004  NA  NA  NA  NA
#> 87  106   0 -0.364963550 -0.364963550 -0.364963550          NA  NA  NA  NA  NA
#> 30  107   0 -1.150449433 -1.150449433 -1.150449433          NA  NA  NA  NA  NA
#> 69  108   0  0.447675407  0.447675407  0.447675407          NA  NA  NA  NA  NA
#> 169 109   1 -0.066396775 -0.066396775 -0.066396775 -0.06639678  NA  NA  NA  NA
#> 202 110   1  0.574880855  0.574880855  0.574880855          NA  NA  NA  NA  NA
#> 196 111   1 -0.391877880 -0.391877880           NA          NA  NA  NA  NA  NA
#> 133 112   1 -0.235937890 -0.235937890 -0.235937890          NA  NA  NA  NA  NA
#> 2   113   0 -0.618695501 -0.618695501 -0.618695501 -0.61869550  NA  NA  NA  NA
#> 186 114   1  0.607833809  0.607833809  0.607833809          NA  NA  NA  NA  NA
#> 208 115   1 -1.930423709 -1.930423709 -1.930423709          NA  NA  NA  NA  NA
#> 35  116   0 -0.172983511 -0.172983511           NA          NA  NA  NA  NA  NA
#> 163 117   1  0.250827790  0.250827790           NA          NA  NA  NA  NA  NA
#> 31  118   0 -1.162498817 -1.162498817           NA          NA  NA  NA  NA  NA
#> 114 119   1 -0.690935395 -0.690935395           NA          NA  NA  NA  NA  NA
#> 57  120   0 -0.234955033 -0.234955033 -0.234955033 -0.23495503  NA  NA  NA  NA
#> 3   121   0 -0.657754469 -0.657754469 -0.657754469          NA  NA  NA  NA  NA
#> 184 122   1 -0.216609744 -0.216609744 -0.216609744          NA  NA  NA  NA  NA
#> 24  123   0 -1.987654126 -1.987654126           NA          NA  NA  NA  NA  NA
#> 76  124   0 -0.354982605 -0.354982605 -0.354982605 -0.35498261  NA  NA  NA  NA
#> 63  125   0  0.848124435  0.848124435  0.848124435          NA  NA  NA  NA  NA
#> 11  126   0  0.868720605  0.868720605  0.868720605          NA  NA  NA  NA  NA
#> 101 127   1 -0.132253303 -0.132253303 -0.132253303          NA  NA  NA  NA  NA
#> 13  128   0 -1.107052754 -1.107052754 -1.107052754          NA  NA  NA  NA  NA
#> 201 129   1  0.472093208  0.472093208           NA          NA  NA  NA  NA  NA
#> 152 130   1 -0.228282333 -0.228282333 -0.228282333          NA  NA  NA  NA  NA
#> 27  131   0 -1.606942139 -1.606942139 -1.606942139          NA  NA  NA  NA  NA
#> 1   132   0 -0.223839052 -0.223839052 -0.223839052          NA  NA  NA  NA  NA
#> 6   133   0  1.239778480  1.239778480  1.239778480  1.23977848  NA  NA  NA  NA
#> 59  134   0  1.190271478  1.190271478  1.190271478  1.19027148  NA  NA  NA  NA
#> 90  135   0  1.375458080  1.375458080  1.375458080          NA  NA  NA  NA  NA
#> 176 136   1  0.086903103  0.086903103  0.086903103  0.08690310  NA  NA  NA  NA
#> 40  137   0  0.367930480  0.367930480           NA          NA  NA  NA  NA  NA
#> 99  138   0  0.182612728  0.182612728           NA          NA  NA  NA  NA  NA
#> 14  139   0  0.391558482  0.391558482  0.391558482          NA  NA  NA  NA  NA
#> 126 140   1  0.748941332  0.748941332           NA          NA  NA  NA  NA  NA
#> 206 141   1 -0.263404679 -0.263404679 -0.263404679          NA  NA  NA  NA  NA
#> 104 142   1  0.402476845  0.402476845  0.402476845          NA  NA  NA  NA  NA
#> 154 143   1 -0.769931534 -0.769931534           NA          NA  NA  NA  NA  NA
#> 151 144   1 -0.430827986 -0.430827986 -0.430827986          NA  NA  NA  NA  NA
#> 139 145   1  0.965517293  0.965517293  0.965517293          NA  NA  NA  NA  NA
#> 91  146   0  2.936759868  2.936759868  2.936759868          NA  NA  NA  NA  NA
#> 36  147   0 -0.539366090 -0.539366090           NA          NA  NA  NA  NA  NA
#> 21  148   0 -0.093649019 -0.093649019           NA          NA  NA  NA  NA  NA
#> 136 149   1  1.291999643  1.291999643           NA          NA  NA  NA  NA  NA
#> 53  150   0  0.483472595  0.483472595           NA          NA  NA  NA  NA  NA
#> 187 151   1  0.000239312  0.000239312           NA          NA  NA  NA  NA  NA
#> 207 152   1 -0.419221617 -0.419221617 -0.419221617 -0.41922162  NA  NA  NA  NA
#> 52  153   0  0.362321932  0.362321932  0.362321932          NA  NA  NA  NA  NA
#> 173 154   1  0.329077694  0.329077694  0.329077694          NA  NA  NA  NA  NA
#> 83  155   0  1.539132115  1.539132115  1.539132115  1.53913211  NA  NA  NA  NA
#> 72  156   0 -0.025271504 -0.025271504           NA          NA  NA  NA  NA  NA
#> 211 157   1 -0.889314955 -0.889314955 -0.889314955          NA  NA  NA  NA  NA
#> 65  158   0 -1.058719062 -1.058719062 -1.058719062          NA  NA  NA  NA  NA
#> 130 159   1  0.002883867  0.002883867  0.002883867          NA  NA  NA  NA  NA
#> 128 160   1 -0.338679322 -0.338679322           NA          NA  NA  NA  NA  NA
#> 178 161   1 -0.013440081 -0.013440081 -0.013440081          NA  NA  NA  NA  NA
#> 145 162   1 -0.605216897 -0.605216897           NA          NA  NA  NA  NA  NA
#> 137 163   1  0.157472272  0.157472272  0.157472272          NA  NA  NA  NA  NA
#> 194 164   1 -0.042452885 -0.042452885           NA          NA  NA  NA  NA  NA
#> 5   165   0 -1.188011877 -1.188011877 -1.188011877          NA  NA  NA  NA  NA
#> 131 166   1  1.283452883  1.283452883  1.283452883  1.28345288  NA  NA  NA  NA
#> 38  167   0  1.024158477  1.024158477  1.024158477          NA  NA  NA  NA  NA
#> 116 168   1  0.532525947  0.532525947  0.532525947          NA  NA  NA  NA  NA
#> 189 169   1 -0.651318395 -0.651318395 -0.651318395          NA  NA  NA  NA  NA
#> 97  170   0  2.732720118  2.732720118  2.732720118          NA  NA  NA  NA  NA
#> 191 171   1  0.190260773  0.190260773           NA          NA  NA  NA  NA  NA
#> 146 172   1  0.495023102  0.495023102           NA          NA  NA  NA  NA  NA
#> 165 173   1  0.291779038  0.291779038  0.291779038          NA  NA  NA  NA  NA
#> 185 174   1 -0.278773118 -0.278773118 -0.278773118 -0.27877312  NA  NA  NA  NA
#> 138 175   1 -0.555270884 -0.555270884 -0.555270884          NA  NA  NA  NA  NA
#> 153 176   1  0.422243636  0.422243636  0.422243636          NA  NA  NA  NA  NA
#> 204 177   1 -1.138775086 -1.138775086 -1.138775086          NA  NA  NA  NA  NA
#> 147 178   1  0.051521241  0.051521241  0.051521241          NA  NA  NA  NA  NA
#> 125 179   1 -0.334662851 -0.334662851 -0.334662851          NA  NA  NA  NA  NA
#> 124 180   1  0.637578012  0.637578012  0.637578012          NA  NA  NA  NA  NA
#> 199 181   1  0.219953958  0.219953958  0.219953958          NA  NA  NA  NA  NA
#> 170 182   1  1.070618491  1.070618491  1.070618491          NA  NA  NA  NA  NA
#> 117 183   1 -0.041284567 -0.041284567           NA          NA  NA  NA  NA  NA
#> 175 184   1  0.194415453  0.194415453  0.194415453          NA  NA  NA  NA  NA
#> 159 185   1  1.982080939  1.982080939  1.982080939          NA  NA  NA  NA  NA
#> 182 186   1 -0.802574794 -0.802574794           NA          NA  NA  NA  NA  NA
#> 20  187   0  0.677665541  0.677665541           NA          NA  NA  NA  NA  NA
#> 181 188   1  0.944153820  0.944153820           NA          NA  NA  NA  NA  NA
#> 109 189   1 -0.461941033 -0.461941033 -0.461941033 -0.46194103  NA  NA  NA  NA
#> 78  190   0  0.768800290  0.768800290           NA          NA  NA  NA  NA  NA
#> 132 191   1  0.193206155  0.193206155  0.193206155          NA  NA  NA  NA  NA
#> 56  192   0  0.306189764  0.306189764  0.306189764          NA  NA  NA  NA  NA
#> 74  193   0  2.224260498  2.224260498           NA          NA  NA  NA  NA  NA
#> 100 194   0  0.262549454  0.262549454  0.262549454          NA  NA  NA  NA  NA
#> 29  195   0  0.594149748  0.594149748           NA          NA  NA  NA  NA  NA
#> 150 196   1  1.065468793  1.065468793  1.065468793          NA  NA  NA  NA  NA
#> 203 197   1  1.707458672  1.707458672  1.707458672          NA  NA  NA  NA  NA
#> 16  198   0  0.294382832  0.294382832  0.294382832  0.29438283  NA  NA  NA  NA
#> 183 199   1 -0.529701796 -0.529701796           NA          NA  NA  NA  NA  NA
#> 85  200   0 -0.236723170 -0.236723170 -0.236723170 -0.23672317  NA  NA  NA  NA
#> 135 201   1  1.034367385  1.034367385  1.034367385          NA  NA  NA  NA  NA
#> 171 202   1  1.043612025  1.043612025           NA          NA  NA  NA  NA  NA
#> 123 203   1  1.330582986  1.330582986           NA          NA  NA  NA  NA  NA
#> 32  204   0  0.432956814  0.432956814           NA          NA  NA  NA  NA  NA
#> 140 205   1 -1.229867169 -1.229867169           NA          NA  NA  NA  NA  NA
#> 112 206   1  0.133645058  0.133645058  0.133645058          NA  NA  NA  NA  NA
#> 148 207   1  3.108110657  3.108110657           NA          NA  NA  NA  NA  NA
#> 198 208   1  0.058599917  0.058599917  0.058599917          NA  NA  NA  NA  NA
#> 180 209   1  2.954902872  2.954902872           NA          NA  NA  NA  NA  NA
#> 58  210   0  1.630736201  1.630736201  1.630736201  1.63073620  NA  NA  NA  NA
#> 177 211   1  0.817882922  0.817882922  0.817882922  0.81788292  NA  NA  NA  NA
#>     X_8 X_9          W_0          W_1          W_2          W_3 W_4 W_5 W_6 W_7
#> 73   NA  NA -0.095487568           NA           NA           NA  NA  NA  NA  NA
#> 89   NA  NA -1.483969090           NA           NA           NA  NA  NA  NA  NA
#> 71   NA  NA -0.984978713           NA           NA           NA  NA  NA  NA  NA
#> 120  NA  NA  0.370679312           NA           NA           NA  NA  NA  NA  NA
#> 77   NA  NA -2.639347932           NA           NA           NA  NA  NA  NA  NA
#> 54   NA  NA -0.866285096           NA           NA           NA  NA  NA  NA  NA
#> 121  NA  NA -0.457168337           NA           NA           NA  NA  NA  NA  NA
#> 141  NA  NA  0.277604643           NA           NA           NA  NA  NA  NA  NA
#> 70   NA  NA -1.454541666           NA           NA           NA  NA  NA  NA  NA
#> 34   NA  NA -0.170030350           NA           NA           NA  NA  NA  NA  NA
#> 50   NA  NA  0.749120907           NA           NA           NA  NA  NA  NA  NA
#> 111  NA  NA -1.390178357           NA           NA           NA  NA  NA  NA  NA
#> 94   NA  NA  1.815126136           NA           NA           NA  NA  NA  NA  NA
#> 107  NA  NA -2.203344553           NA           NA           NA  NA  NA  NA  NA
#> 19   NA  NA -0.099402677           NA           NA           NA  NA  NA  NA  NA
#> 115  NA  NA  0.245524670           NA           NA           NA  NA  NA  NA  NA
#> 23   NA  NA  0.217880795 -0.145100465           NA           NA  NA  NA  NA  NA
#> 47   NA  NA -1.075999051           NA           NA           NA  NA  NA  NA  NA
#> 67   NA  NA  0.002995459           NA           NA           NA  NA  NA  NA  NA
#> 81   NA  NA -1.231787034           NA           NA           NA  NA  NA  NA  NA
#> 142  NA  NA  0.212185450           NA           NA           NA  NA  NA  NA  NA
#> 37   NA  NA  1.122275564           NA           NA           NA  NA  NA  NA  NA
#> 8    NA  NA -0.397932510           NA           NA           NA  NA  NA  NA  NA
#> 79   NA  NA  0.181176517 -0.634436761           NA           NA  NA  NA  NA  NA
#> 95   NA  NA -0.743724035 -0.428348752           NA           NA  NA  NA  NA  NA
#> 26   NA  NA  0.460046529           NA           NA           NA  NA  NA  NA  NA
#> 155  NA  NA  1.518396696  1.994615842           NA           NA  NA  NA  NA  NA
#> 48   NA  NA  0.712021433           NA           NA           NA  NA  NA  NA  NA
#> 158  NA  NA -2.672681664           NA           NA           NA  NA  NA  NA  NA
#> 200  NA  NA -0.623431035           NA           NA           NA  NA  NA  NA  NA
#> 134  NA  NA  0.278552309           NA           NA           NA  NA  NA  NA  NA
#> 156  NA  NA  0.004725478           NA           NA           NA  NA  NA  NA  NA
#> 119  NA  NA -0.193716472           NA           NA           NA  NA  NA  NA  NA
#> 127  NA  NA  0.331400128 -0.953401170           NA           NA  NA  NA  NA  NA
#> 103  NA  NA  0.026496136           NA           NA           NA  NA  NA  NA  NA
#> 193  NA  NA -0.700733592 -0.966037492           NA           NA  NA  NA  NA  NA
#> 118  NA  NA -0.399372778           NA           NA           NA  NA  NA  NA  NA
#> 209  NA  NA  0.443082681           NA           NA           NA  NA  NA  NA  NA
#> 98   NA  NA -0.363756203           NA           NA           NA  NA  NA  NA  NA
#> 61   NA  NA -1.149284301           NA           NA           NA  NA  NA  NA  NA
#> 4    NA  NA -1.096098481           NA           NA           NA  NA  NA  NA  NA
#> 157  NA  NA -2.197902005           NA           NA           NA  NA  NA  NA  NA
#> 7    NA  NA -1.156774088           NA           NA           NA  NA  NA  NA  NA
#> 66   NA  NA -0.649366033           NA           NA           NA  NA  NA  NA  NA
#> 60   NA  NA  1.193529373  0.190199197           NA           NA  NA  NA  NA  NA
#> 43   NA  NA -2.325448321 -1.777459761           NA           NA  NA  NA  NA  NA
#> 210  NA  NA -1.682346208           NA           NA           NA  NA  NA  NA  NA
#> 28   NA  NA  0.229830034 -0.390359715           NA           NA  NA  NA  NA  NA
#> 17   NA  NA  0.114119709  0.959147417           NA           NA  NA  NA  NA  NA
#> 162  NA  NA -1.762610526           NA           NA           NA  NA  NA  NA  NA
#> 205  NA  NA  0.886664253           NA           NA           NA  NA  NA  NA  NA
#> 80   NA  NA -0.659518318           NA           NA           NA  NA  NA  NA  NA
#> 188  NA  NA  0.706512090  0.112258388           NA           NA  NA  NA  NA  NA
#> 25   NA  NA -0.729974507           NA           NA           NA  NA  NA  NA  NA
#> 84   NA  NA  0.512433788  0.797193998           NA           NA  NA  NA  NA  NA
#> 108  NA  NA  0.256595324  0.330112784 -0.498800984           NA  NA  NA  NA  NA
#> 9    NA  NA -1.606157999           NA           NA           NA  NA  NA  NA  NA
#> 144  NA  NA -2.153929733           NA           NA           NA  NA  NA  NA  NA
#> 195  NA  NA -0.066677906  0.717702714           NA           NA  NA  NA  NA  NA
#> 168  NA  NA  0.289605229  0.542413481           NA           NA  NA  NA  NA  NA
#> 149  NA  NA  0.228321508  0.891190434 -1.316155745           NA  NA  NA  NA  NA
#> 41   NA  NA  1.029129162 -0.464424418           NA           NA  NA  NA  NA  NA
#> 113  NA  NA  0.534821963  1.181781255           NA           NA  NA  NA  NA  NA
#> 167  NA  NA  0.358871420  0.700483815 -0.564814549           NA  NA  NA  NA  NA
#> 44   NA  NA  0.667368499 -0.839747508           NA           NA  NA  NA  NA  NA
#> 166  NA  NA -1.665574147 -0.363812033           NA           NA  NA  NA  NA  NA
#> 55   NA  NA -0.472669441 -1.017509154           NA           NA  NA  NA  NA  NA
#> 92   NA  NA -0.442105717 -0.103325191           NA           NA  NA  NA  NA  NA
#> 197  NA  NA  1.136053954  1.685179557  0.553974501           NA  NA  NA  NA  NA
#> 143  NA  NA  2.074215381  1.251233285  0.777943247           NA  NA  NA  NA  NA
#> 12   NA  NA -0.017609808 -1.114157044           NA           NA  NA  NA  NA  NA
#> 164  NA  NA  1.981306546  1.166318665  0.151650132           NA  NA  NA  NA  NA
#> 22   NA  NA -0.955102753 -0.782831639           NA           NA  NA  NA  NA  NA
#> 42   NA  NA  0.231004500 -1.553093354 -1.226790624           NA  NA  NA  NA  NA
#> 46   NA  NA  0.939843814  0.808411757           NA           NA  NA  NA  NA  NA
#> 68   NA  NA -0.722859130 -1.510160948           NA           NA  NA  NA  NA  NA
#> 18   NA  NA  0.214645522 -0.069511408           NA           NA  NA  NA  NA  NA
#> 179  NA  NA  2.176927062  1.186101557  1.151982609           NA  NA  NA  NA  NA
#> 75   NA  NA  1.140454283  0.830653820  0.664488832           NA  NA  NA  NA  NA
#> 96   NA  NA -0.606809943 -0.845820087 -1.306204414 -1.497200025  NA  NA  NA  NA
#> 49   NA  NA  0.039857822  0.105129280           NA           NA  NA  NA  NA  NA
#> 33   NA  NA  0.214883569 -0.506676632 -0.884221563           NA  NA  NA  NA  NA
#> 129  NA  NA  1.000724557  0.765147460           NA           NA  NA  NA  NA  NA
#> 122  NA  NA -2.603336740 -2.188518422           NA           NA  NA  NA  NA  NA
#> 102  NA  NA  0.698129860  0.596557627  0.854427448 -0.006484826  NA  NA  NA  NA
#> 51   NA  NA -0.038730879  0.272490237  1.218097915  0.032757062  NA  NA  NA  NA
#> 88   NA  NA -1.425318914 -0.657828091           NA           NA  NA  NA  NA  NA
#> 161  NA  NA  1.191279635  0.493914632           NA           NA  NA  NA  NA  NA
#> 160  NA  NA -0.513738815  1.029047751 -0.122231761           NA  NA  NA  NA  NA
#> 105  NA  NA -0.792793107 -1.507834756           NA           NA  NA  NA  NA  NA
#> 86   NA  NA -0.247629059 -0.885155679           NA           NA  NA  NA  NA  NA
#> 10   NA  NA  0.469382921 -0.333735481 -0.216058216  0.188390500  NA  NA  NA  NA
#> 15   NA  NA  0.527380537  1.045898034  1.406757767           NA  NA  NA  NA  NA
#> 93   NA  NA  0.511908212  0.803981424  1.078027280 -2.093849242  NA  NA  NA  NA
#> 64   NA  NA  1.782626476  0.327536272 -0.764155521           NA  NA  NA  NA  NA
#> 110  NA  NA -0.009079277  1.709899033  0.300845824  1.141681956  NA  NA  NA  NA
#> 62   NA  NA  1.026890967 -0.005846515           NA           NA  NA  NA  NA  NA
#> 172  NA  NA -0.014161976 -0.123027718  0.570356577 -0.444381092  NA  NA  NA  NA
#> 174  NA  NA  0.847661951  0.183682963 -1.250677267           NA  NA  NA  NA  NA
#> 192  NA  NA  0.763149796 -0.877896333  0.590447192           NA  NA  NA  NA  NA
#> 106  NA  NA  0.478843397 -0.329900101  0.675034534           NA  NA  NA  NA  NA
#> 82   NA  NA  0.065509091 -1.739229057 -1.764678685 -0.201182410  NA  NA  NA  NA
#> 190  NA  NA -0.842275576 -0.210158144  1.048641320  0.424402859  NA  NA  NA  NA
#> 39   NA  NA  1.617593486  1.572118110  1.461140687  0.519726167  NA  NA  NA  NA
#> 45   NA  NA -0.933266617 -1.300056442 -0.969843889 -0.409563863  NA  NA  NA  NA
#> 87   NA  NA  0.080245600  0.593535586  0.069786523           NA  NA  NA  NA  NA
#> 30   NA  NA  0.291347371  1.452500829  1.054807663           NA  NA  NA  NA  NA
#> 69   NA  NA -1.903557577 -0.965697542 -1.599225756           NA  NA  NA  NA  NA
#> 169  NA  NA -0.904853451 -0.827886530 -0.375117579 -1.722072868  NA  NA  NA  NA
#> 202  NA  NA -0.350173032 -1.023820959  0.857207961           NA  NA  NA  NA  NA
#> 196  NA  NA  0.926434516  0.956423145           NA           NA  NA  NA  NA  NA
#> 133  NA  NA  1.190973428 -0.373016726 -0.772648365           NA  NA  NA  NA  NA
#> 2    NA  NA  1.056669487 -0.772981776 -0.513218493 -0.002046591  NA  NA  NA  NA
#> 186  NA  NA  0.705479091  0.153873055 -0.965838602           NA  NA  NA  NA  NA
#> 208  NA  NA  0.828217986  0.142049468 -1.397427088           NA  NA  NA  NA  NA
#> 35   NA  NA -1.146577802  0.983944967           NA           NA  NA  NA  NA  NA
#> 163  NA  NA  0.070864780  0.647170694           NA           NA  NA  NA  NA  NA
#> 31   NA  NA -1.041109234  0.897545024           NA           NA  NA  NA  NA  NA
#> 114  NA  NA -0.853781551  0.045809490           NA           NA  NA  NA  NA  NA
#> 57   NA  NA  2.010537300  0.765827347  0.513090275  0.415213847  NA  NA  NA  NA
#> 3    NA  NA -0.818070273  0.257131900 -0.146771542           NA  NA  NA  NA  NA
#> 184  NA  NA -1.179172756 -0.794648100  0.006680216           NA  NA  NA  NA  NA
#> 24   NA  NA  0.870014181  0.361632450           NA           NA  NA  NA  NA  NA
#> 76   NA  NA  0.405993998  0.346493533  1.290332568  0.595019085  NA  NA  NA  NA
#> 63   NA  NA  0.129228299  0.122788989 -0.783605787           NA  NA  NA  NA  NA
#> 11   NA  NA -1.340219960 -1.306432477 -1.171477421           NA  NA  NA  NA  NA
#> 101  NA  NA  0.851248488 -0.960536609  0.179833702           NA  NA  NA  NA  NA
#> 13   NA  NA -0.112343567 -0.160108356 -1.200064483           NA  NA  NA  NA  NA
#> 201  NA  NA -2.315046633 -0.396754751           NA           NA  NA  NA  NA  NA
#> 152  NA  NA  1.630881040  0.806971057  0.657808273           NA  NA  NA  NA  NA
#> 27   NA  NA  1.620908680  1.012545650  0.423242191           NA  NA  NA  NA  NA
#> 1    NA  NA  1.292812744  1.479563989  1.635903604           NA  NA  NA  NA  NA
#> 6    NA  NA  1.469164347 -0.198018407 -0.452565743 -0.015797535  NA  NA  NA  NA
#> 59   NA  NA  1.205881072 -0.115210171  0.118662628  1.903628559  NA  NA  NA  NA
#> 90   NA  NA -0.730853747 -1.046695306 -1.082141583           NA  NA  NA  NA  NA
#> 176  NA  NA -0.353964586 -1.411409803 -0.954897336 -1.689800190  NA  NA  NA  NA
#> 40   NA  NA -0.449310607 -0.186063990           NA           NA  NA  NA  NA  NA
#> 99   NA  NA -0.540114035  0.890428658           NA           NA  NA  NA  NA  NA
#> 14   NA  NA -0.886896115  0.505491655 -0.058148876           NA  NA  NA  NA  NA
#> 126  NA  NA -0.047750442  0.470548422           NA           NA  NA  NA  NA  NA
#> 206  NA  NA -0.866087336  0.239593708  0.133353685           NA  NA  NA  NA  NA
#> 104  NA  NA -1.013155029 -1.707372777  0.378516311           NA  NA  NA  NA  NA
#> 154  NA  NA  0.574902490  0.484564487           NA           NA  NA  NA  NA  NA
#> 151  NA  NA  1.340991631  0.793857097 -0.038277348           NA  NA  NA  NA  NA
#> 139  NA  NA  0.076758322  0.344810450 -0.191108988           NA  NA  NA  NA  NA
#> 91   NA  NA -1.854374754 -1.024557857 -0.876148543           NA  NA  NA  NA  NA
#> 36   NA  NA -0.941360999  1.027955780           NA           NA  NA  NA  NA  NA
#> 21   NA  NA -0.026483949  0.107372181           NA           NA  NA  NA  NA  NA
#> 136  NA  NA  0.787901015  0.726023557           NA           NA  NA  NA  NA  NA
#> 53   NA  NA  0.228799037  0.092340630           NA           NA  NA  NA  NA  NA
#> 187  NA  NA -0.754839909  0.456029548           NA           NA  NA  NA  NA  NA
#> 207  NA  NA  0.042614733 -0.201076025 -0.174584900  0.178795103  NA  NA  NA  NA
#> 52   NA  NA  0.286257000  0.836470166  1.269907956           NA  NA  NA  NA  NA
#> 173  NA  NA -1.424230301  0.798829241 -0.605277583           NA  NA  NA  NA  NA
#> 83   NA  NA -1.317665756 -1.112540541 -0.386540554  0.331799644  NA  NA  NA  NA
#> 72   NA  NA -0.890257399  0.820049130           NA           NA  NA  NA  NA  NA
#> 211  NA  NA  1.085570464  2.325220498  1.884016488           NA  NA  NA  NA  NA
#> 65   NA  NA  0.086101145  1.449681390  0.815446208           NA  NA  NA  NA  NA
#> 130  NA  NA  0.293359195  1.050134106  0.454050912           NA  NA  NA  NA  NA
#> 128  NA  NA -0.272712742  0.742394690           NA           NA  NA  NA  NA  NA
#> 178  NA  NA  0.520021399  0.114621910  1.271470219           NA  NA  NA  NA  NA
#> 145  NA  NA  0.216480811  0.108924734           NA           NA  NA  NA  NA  NA
#> 137  NA  NA  0.576398201  0.429029574  1.244996358           NA  NA  NA  NA  NA
#> 194  NA  NA -0.394602879 -0.642092390           NA           NA  NA  NA  NA  NA
#> 5    NA  NA  0.507616069  0.495481465  1.561817465           NA  NA  NA  NA  NA
#> 131  NA  NA -0.234121872  1.769773932  1.834072200  0.827317973  NA  NA  NA  NA
#> 38   NA  NA  0.221598954  1.894852489  0.920546211           NA  NA  NA  NA  NA
#> 116  NA  NA -1.129847722 -1.224026051  0.136613075           NA  NA  NA  NA  NA
#> 189  NA  NA  1.910394287  1.791213228  1.676404890           NA  NA  NA  NA  NA
#> 97   NA  NA -0.108392294 -0.990890694 -1.764973150           NA  NA  NA  NA  NA
#> 191  NA  NA -0.159600759  1.267249887           NA           NA  NA  NA  NA  NA
#> 146  NA  NA  1.497935962  1.798082307           NA           NA  NA  NA  NA  NA
#> 165  NA  NA  0.548449812  0.110769206  0.981880866           NA  NA  NA  NA  NA
#> 185  NA  NA  0.047049887  1.086827859  0.241362405 -0.465850161  NA  NA  NA  NA
#> 138  NA  NA -0.251747207  1.327477435  1.013379159           NA  NA  NA  NA  NA
#> 153  NA  NA  1.002425646  1.450728387  1.732756207           NA  NA  NA  NA  NA
#> 204  NA  NA  2.107249939  1.769857490  2.652235961           NA  NA  NA  NA  NA
#> 147  NA  NA -0.348442242  0.794341356 -0.469047053           NA  NA  NA  NA  NA
#> 125  NA  NA  0.949086618  0.552284267  0.374995053           NA  NA  NA  NA  NA
#> 124  NA  NA  1.648359741  2.242678600  2.162234073           NA  NA  NA  NA  NA
#> 199  NA  NA -2.303321618 -0.081368703  0.190855510           NA  NA  NA  NA  NA
#> 170  NA  NA  0.476040369 -0.544946055  0.539172256           NA  NA  NA  NA  NA
#> 117  NA  NA  0.108344013  0.218656992           NA           NA  NA  NA  NA  NA
#> 175  NA  NA  0.054490877  0.752823809  1.119287739           NA  NA  NA  NA  NA
#> 159  NA  NA  0.931079174  0.553794579 -0.911936085           NA  NA  NA  NA  NA
#> 182  NA  NA  0.516508078  1.298886574           NA           NA  NA  NA  NA  NA
#> 20   NA  NA  0.517502546 -0.137006008           NA           NA  NA  NA  NA  NA
#> 181  NA  NA  0.075108124 -0.603231137           NA           NA  NA  NA  NA  NA
#> 109  NA  NA  0.412373712  0.496681871  0.255663956  1.731879596  NA  NA  NA  NA
#> 78   NA  NA  1.421992573  1.542504675           NA           NA  NA  NA  NA  NA
#> 132  NA  NA -0.573606808  1.023881489  2.201360543           NA  NA  NA  NA  NA
#> 56   NA  NA -0.021299324 -0.626894412  0.947737591           NA  NA  NA  NA  NA
#> 74   NA  NA -1.558264712  0.523298078           NA           NA  NA  NA  NA  NA
#> 100  NA  NA -0.534797772 -0.009829929  1.198475953           NA  NA  NA  NA  NA
#> 29   NA  NA -1.154317933 -0.879938872           NA           NA  NA  NA  NA  NA
#> 150  NA  NA  0.849416980  0.828604184  1.612586109           NA  NA  NA  NA  NA
#> 203  NA  NA  1.173562215  1.108747313  2.295233105           NA  NA  NA  NA  NA
#> 16   NA  NA -0.118895848  0.206800373 -0.342692399  0.804242037  NA  NA  NA  NA
#> 183  NA  NA -0.425578710 -0.803443811           NA           NA  NA  NA  NA  NA
#> 85   NA  NA -0.481200350  2.362095485  0.980401806  0.950186249  NA  NA  NA  NA
#> 135  NA  NA -1.229962916 -0.946372304 -1.440901213           NA  NA  NA  NA  NA
#> 171  NA  NA -0.384014961 -0.555775125           NA           NA  NA  NA  NA  NA
#> 123  NA  NA  0.516811277 -0.272089290           NA           NA  NA  NA  NA  NA
#> 32   NA  NA  1.209559056  0.457896323           NA           NA  NA  NA  NA  NA
#> 140  NA  NA -0.472478384 -0.522151935           NA           NA  NA  NA  NA  NA
#> 112  NA  NA -1.249604658 -0.543871033  0.291213408           NA  NA  NA  NA  NA
#> 148  NA  NA  1.264060940 -0.571821422           NA           NA  NA  NA  NA  NA
#> 198  NA  NA  0.102162280  1.098382294 -0.574262680           NA  NA  NA  NA  NA
#> 180  NA  NA  1.599553817  0.575481861           NA           NA  NA  NA  NA  NA
#> 58   NA  NA  0.860095869  0.410377915  1.558343492  1.175705030  NA  NA  NA  NA
#> 177  NA  NA -0.573861712  0.134849957 -0.570150621 -0.739316247  NA  NA  NA  NA
#>     W_8 W_9          L_0          L_1          L_2         L_3 L_4 L_5 L_6 L_7
#> 73   NA  NA -0.794908725           NA           NA          NA  NA  NA  NA  NA
#> 89   NA  NA  1.648630293           NA           NA          NA  NA  NA  NA  NA
#> 71   NA  NA  0.009935887           NA           NA          NA  NA  NA  NA  NA
#> 120  NA  NA  1.525530898           NA           NA          NA  NA  NA  NA  NA
#> 77   NA  NA  0.010827752           NA           NA          NA  NA  NA  NA  NA
#> 54   NA  NA -1.196253657           NA           NA          NA  NA  NA  NA  NA
#> 121  NA  NA -1.115205699           NA           NA          NA  NA  NA  NA  NA
#> 141  NA  NA -0.556643307           NA           NA          NA  NA  NA  NA  NA
#> 70   NA  NA  1.104915317           NA           NA          NA  NA  NA  NA  NA
#> 34   NA  NA  0.503835290           NA           NA          NA  NA  NA  NA  NA
#> 50   NA  NA -1.197744279           NA           NA          NA  NA  NA  NA  NA
#> 111  NA  NA  1.058795806           NA           NA          NA  NA  NA  NA  NA
#> 94   NA  NA  1.164476460           NA           NA          NA  NA  NA  NA  NA
#> 107  NA  NA  1.480581378           NA           NA          NA  NA  NA  NA  NA
#> 19   NA  NA  0.651809529           NA           NA          NA  NA  NA  NA  NA
#> 115  NA  NA  0.677606918           NA           NA          NA  NA  NA  NA  NA
#> 23   NA  NA -0.382772108  1.847165235           NA          NA  NA  NA  NA  NA
#> 47   NA  NA -0.314709642           NA           NA          NA  NA  NA  NA  NA
#> 67   NA  NA -0.313104798           NA           NA          NA  NA  NA  NA  NA
#> 81   NA  NA  2.094910336           NA           NA          NA  NA  NA  NA  NA
#> 142  NA  NA  0.956033786           NA           NA          NA  NA  NA  NA  NA
#> 37   NA  NA -2.123720123           NA           NA          NA  NA  NA  NA  NA
#> 8    NA  NA  1.802232725           NA           NA          NA  NA  NA  NA  NA
#> 79   NA  NA -0.278958962  0.156648956           NA          NA  NA  NA  NA  NA
#> 95   NA  NA  0.378352284 -1.534038035           NA          NA  NA  NA  NA  NA
#> 26   NA  NA -0.969804695           NA           NA          NA  NA  NA  NA  NA
#> 155  NA  NA  0.637275130  0.023343978           NA          NA  NA  NA  NA  NA
#> 48   NA  NA -2.270163564           NA           NA          NA  NA  NA  NA  NA
#> 158  NA  NA -1.477314055           NA           NA          NA  NA  NA  NA  NA
#> 200  NA  NA -1.155078015           NA           NA          NA  NA  NA  NA  NA
#> 134  NA  NA -0.760119705           NA           NA          NA  NA  NA  NA  NA
#> 156  NA  NA -0.855751849           NA           NA          NA  NA  NA  NA  NA
#> 119  NA  NA -0.194082531           NA           NA          NA  NA  NA  NA  NA
#> 127  NA  NA  1.441918671  1.030138114           NA          NA  NA  NA  NA  NA
#> 103  NA  NA  1.102330933           NA           NA          NA  NA  NA  NA  NA
#> 193  NA  NA -0.874565020 -0.433380669           NA          NA  NA  NA  NA  NA
#> 118  NA  NA -1.821175364           NA           NA          NA  NA  NA  NA  NA
#> 209  NA  NA  0.914832753           NA           NA          NA  NA  NA  NA  NA
#> 98   NA  NA -0.533248173           NA           NA          NA  NA  NA  NA  NA
#> 61   NA  NA  2.128569801           NA           NA          NA  NA  NA  NA  NA
#> 4    NA  NA  2.032670306           NA           NA          NA  NA  NA  NA  NA
#> 157  NA  NA  0.177990132           NA           NA          NA  NA  NA  NA  NA
#> 7    NA  NA  1.141601107           NA           NA          NA  NA  NA  NA  NA
#> 66   NA  NA  0.098571479           NA           NA          NA  NA  NA  NA  NA
#> 60   NA  NA -1.132853826 -1.235585284           NA          NA  NA  NA  NA  NA
#> 43   NA  NA -0.053452190  0.512907378           NA          NA  NA  NA  NA  NA
#> 210  NA  NA  0.438172274           NA           NA          NA  NA  NA  NA  NA
#> 28   NA  NA  1.497799918  0.854217903           NA          NA  NA  NA  NA  NA
#> 17   NA  NA  0.222469688  1.383109093           NA          NA  NA  NA  NA  NA
#> 162  NA  NA  0.645166675           NA           NA          NA  NA  NA  NA  NA
#> 205  NA  NA  1.462052583           NA           NA          NA  NA  NA  NA  NA
#> 80   NA  NA -1.498630983           NA           NA          NA  NA  NA  NA  NA
#> 188  NA  NA  0.740403519  1.274302001           NA          NA  NA  NA  NA  NA
#> 25   NA  NA  0.020275989           NA           NA          NA  NA  NA  NA  NA
#> 84   NA  NA -0.659231022 -0.507741402           NA          NA  NA  NA  NA  NA
#> 108  NA  NA  0.571464822  0.720579915  0.039870688          NA  NA  NA  NA  NA
#> 9    NA  NA -1.759947009           NA           NA          NA  NA  NA  NA  NA
#> 144  NA  NA  0.902280636           NA           NA          NA  NA  NA  NA  NA
#> 195  NA  NA  0.495780891  0.211410102           NA          NA  NA  NA  NA  NA
#> 168  NA  NA -0.763095428 -0.189107011           NA          NA  NA  NA  NA  NA
#> 149  NA  NA  0.111138463  0.389695077  0.263024244          NA  NA  NA  NA  NA
#> 41   NA  NA  0.512247411  1.619882170           NA          NA  NA  NA  NA  NA
#> 113  NA  NA -0.325855498  0.535951243           NA          NA  NA  NA  NA  NA
#> 167  NA  NA  0.689875751  0.703713479  0.426788292          NA  NA  NA  NA  NA
#> 44   NA  NA -1.126027190 -1.432368367           NA          NA  NA  NA  NA  NA
#> 166  NA  NA -0.464071382  0.627261293           NA          NA  NA  NA  NA  NA
#> 55   NA  NA  0.153955408 -0.101929018           NA          NA  NA  NA  NA  NA
#> 92   NA  NA -0.124810495 -0.007107784           NA          NA  NA  NA  NA  NA
#> 197  NA  NA  0.216581911 -1.190064460  0.490368922          NA  NA  NA  NA  NA
#> 143  NA  NA  1.297360843  1.116642766 -0.148035420          NA  NA  NA  NA  NA
#> 12   NA  NA  0.180525062 -0.252004509           NA          NA  NA  NA  NA  NA
#> 164  NA  NA -1.849473355 -1.672560384 -0.753162179          NA  NA  NA  NA  NA
#> 22   NA  NA -0.145914654  1.687191684           NA          NA  NA  NA  NA  NA
#> 42   NA  NA -0.524668998  1.090507677  0.438490393          NA  NA  NA  NA  NA
#> 46   NA  NA  1.456810877  1.450714957           NA          NA  NA  NA  NA  NA
#> 68   NA  NA  1.318761959 -0.109363251           NA          NA  NA  NA  NA  NA
#> 18   NA  NA  0.386555905  1.224726357           NA          NA  NA  NA  NA  NA
#> 179  NA  NA -0.886857879  0.133626582  0.927508093          NA  NA  NA  NA  NA
#> 75   NA  NA -0.771209291 -1.540592499 -0.929909731          NA  NA  NA  NA  NA
#> 96   NA  NA -0.562968455 -0.631492027 -0.635850905 -0.86242486  NA  NA  NA  NA
#> 49   NA  NA -0.072309982 -0.925829711           NA          NA  NA  NA  NA  NA
#> 33   NA  NA  0.116506419 -1.437663340 -0.883531064          NA  NA  NA  NA  NA
#> 129  NA  NA  0.123690580  0.415690684           NA          NA  NA  NA  NA  NA
#> 122  NA  NA  1.014189038  0.046545844           NA          NA  NA  NA  NA  NA
#> 102  NA  NA  1.537533686  0.636537766  0.402181441  0.78957469  NA  NA  NA  NA
#> 51   NA  NA -1.671914220 -1.822887636 -1.148980889 -1.23820212  NA  NA  NA  NA
#> 88   NA  NA -1.317752176 -0.222507614           NA          NA  NA  NA  NA  NA
#> 161  NA  NA  2.160654795  2.320931135           NA          NA  NA  NA  NA  NA
#> 160  NA  NA  0.983978841 -0.314836683 -0.601188181          NA  NA  NA  NA  NA
#> 105  NA  NA -0.272388414  0.568849890           NA          NA  NA  NA  NA  NA
#> 86   NA  NA -0.222495148  0.948454768           NA          NA  NA  NA  NA  NA
#> 10   NA  NA -0.147538958  0.905021833 -0.223729473  0.03208246  NA  NA  NA  NA
#> 15   NA  NA -0.017123892  0.298582762  1.062457222          NA  NA  NA  NA  NA
#> 93   NA  NA  0.843557918 -0.721936314 -0.755508764 -0.84319229  NA  NA  NA  NA
#> 64   NA  NA  0.420789951 -0.475192101 -0.840682068          NA  NA  NA  NA  NA
#> 110  NA  NA -1.652279469  0.463488203 -1.180998537 -1.24761156  NA  NA  NA  NA
#> 62   NA  NA -0.633371120  0.518692063           NA          NA  NA  NA  NA  NA
#> 172  NA  NA -0.276593368 -1.031329602  0.494473021  0.32852040  NA  NA  NA  NA
#> 174  NA  NA -1.302929075 -0.393432494 -1.063808354          NA  NA  NA  NA  NA
#> 192  NA  NA -1.526714454 -0.398992354 -0.891171381          NA  NA  NA  NA  NA
#> 106  NA  NA -0.592491787  0.360147587 -0.411436240          NA  NA  NA  NA  NA
#> 82   NA  NA -1.266096352 -0.408157957  1.028442068  0.03108140  NA  NA  NA  NA
#> 190  NA  NA -0.827964487 -0.352128637 -1.095101121 -1.35013121  NA  NA  NA  NA
#> 39   NA  NA  0.446124915  1.148293433  0.380052177 -0.24630585  NA  NA  NA  NA
#> 45   NA  NA  0.066932290 -1.950752636  0.802160564  0.33511934  NA  NA  NA  NA
#> 87   NA  NA  0.761639720  0.219456755  0.352958834          NA  NA  NA  NA  NA
#> 30   NA  NA  0.944487137  0.053714977 -0.071169966          NA  NA  NA  NA  NA
#> 69   NA  NA -1.956659025 -1.406519220 -1.734844053          NA  NA  NA  NA  NA
#> 169  NA  NA -0.519461812  1.138382791  1.072933586  0.40506514  NA  NA  NA  NA
#> 202  NA  NA  2.028201820  1.127551055  1.758392538          NA  NA  NA  NA  NA
#> 196  NA  NA  0.980772093  1.573857674           NA          NA  NA  NA  NA  NA
#> 133  NA  NA  1.475938321  1.381269804  2.257525416          NA  NA  NA  NA  NA
#> 2    NA  NA -0.159859382  0.646887889  0.025420368 -2.09507253  NA  NA  NA  NA
#> 186  NA  NA -0.862364263 -0.589497305 -0.582412077          NA  NA  NA  NA  NA
#> 208  NA  NA  0.228503882 -0.493365053 -0.217950695          NA  NA  NA  NA  NA
#> 35   NA  NA -0.760092491 -1.805426437           NA          NA  NA  NA  NA  NA
#> 163  NA  NA -1.581449157 -1.636674279           NA          NA  NA  NA  NA  NA
#> 31   NA  NA  1.064098669  0.770281653           NA          NA  NA  NA  NA  NA
#> 114  NA  NA  0.452656967  0.126801466           NA          NA  NA  NA  NA  NA
#> 57   NA  NA -0.966658739 -1.857603202 -1.625695619 -2.28943977  NA  NA  NA  NA
#> 3    NA  NA  0.782547913 -1.470561444  1.867836783          NA  NA  NA  NA  NA
#> 184  NA  NA -0.736173249  0.079815112 -0.716390964          NA  NA  NA  NA  NA
#> 24   NA  NA  0.501218212  1.290334441           NA          NA  NA  NA  NA  NA
#> 76   NA  NA  2.060636237  2.224246055  3.117867271  1.73071096  NA  NA  NA  NA
#> 63   NA  NA  0.605264141 -1.029007487 -1.802418921          NA  NA  NA  NA  NA
#> 11   NA  NA  0.917805322  0.244589874 -1.202573110          NA  NA  NA  NA  NA
#> 101  NA  NA -0.381741142 -0.788186265  1.500498837          NA  NA  NA  NA  NA
#> 13   NA  NA -1.638977166 -0.531805789  0.825987244          NA  NA  NA  NA  NA
#> 201  NA  NA -2.182697279 -1.222826894           NA          NA  NA  NA  NA  NA
#> 152  NA  NA -0.560786061 -1.879958321  0.139604637          NA  NA  NA  NA  NA
#> 27   NA  NA -0.084679412 -0.619092719 -0.740293696          NA  NA  NA  NA  NA
#> 1    NA  NA -1.891566677 -0.804077972 -1.136570274          NA  NA  NA  NA  NA
#> 6    NA  NA  0.084687426  0.869638414 -0.007844704  1.41924980  NA  NA  NA  NA
#> 59   NA  NA -2.185488606 -1.619178570 -0.914789429 -1.01517082  NA  NA  NA  NA
#> 90   NA  NA -0.301904177  0.374152164  0.635813852          NA  NA  NA  NA  NA
#> 176  NA  NA -2.820126313 -1.342183327 -1.296127137 -1.64553382  NA  NA  NA  NA
#> 40   NA  NA  0.086083269  0.319532233           NA          NA  NA  NA  NA  NA
#> 99   NA  NA -0.046764494 -0.853466898           NA          NA  NA  NA  NA  NA
#> 14   NA  NA  1.397524883  1.008869754  1.041511067          NA  NA  NA  NA  NA
#> 126  NA  NA  1.835764370 -0.242203642           NA          NA  NA  NA  NA  NA
#> 206  NA  NA  1.662589287  1.423176525  1.066804857          NA  NA  NA  NA  NA
#> 104  NA  NA  0.745789066 -0.439370457 -0.464676965          NA  NA  NA  NA  NA
#> 154  NA  NA -0.740002546  0.341021898           NA          NA  NA  NA  NA  NA
#> 151  NA  NA -0.997307241 -1.310502719 -1.242959721          NA  NA  NA  NA  NA
#> 139  NA  NA  0.375745378 -0.771261866 -0.283518832          NA  NA  NA  NA  NA
#> 91   NA  NA -1.196621759 -0.539973827 -1.456259195          NA  NA  NA  NA  NA
#> 36   NA  NA  1.915441266  1.282980487           NA          NA  NA  NA  NA  NA
#> 21   NA  NA -0.388766460 -0.145335443           NA          NA  NA  NA  NA  NA
#> 136  NA  NA  1.380071385  0.080080154           NA          NA  NA  NA  NA  NA
#> 53   NA  NA -0.378212300 -1.069585832           NA          NA  NA  NA  NA  NA
#> 187  NA  NA -0.954110761 -0.115713770           NA          NA  NA  NA  NA  NA
#> 207  NA  NA  0.385244932  0.492436120  0.327263524  0.95025561  NA  NA  NA  NA
#> 52   NA  NA -0.327353734  0.339400827  0.674134008          NA  NA  NA  NA  NA
#> 173  NA  NA -0.323890827 -1.470802305 -0.178335109          NA  NA  NA  NA  NA
#> 83   NA  NA  0.759673191  0.593985266  0.364475519  0.61401939  NA  NA  NA  NA
#> 72   NA  NA  0.945645245  1.297472890           NA          NA  NA  NA  NA  NA
#> 211  NA  NA -0.359191132 -1.367071861  0.098374159          NA  NA  NA  NA  NA
#> 65   NA  NA -2.120953225 -0.953874493 -1.831758661          NA  NA  NA  NA  NA
#> 130  NA  NA -1.066299400 -1.531253617 -1.509590545          NA  NA  NA  NA  NA
#> 128  NA  NA  2.385231987  2.993528520           NA          NA  NA  NA  NA  NA
#> 178  NA  NA  0.218882485 -0.153976120  0.871052396          NA  NA  NA  NA  NA
#> 145  NA  NA -0.784392191 -0.321718281           NA          NA  NA  NA  NA  NA
#> 137  NA  NA  0.771410213  1.424297179  0.447056608          NA  NA  NA  NA  NA
#> 194  NA  NA -1.085479097 -1.531818574           NA          NA  NA  NA  NA  NA
#> 5    NA  NA  0.892913949  1.381916805 -0.019146638          NA  NA  NA  NA  NA
#> 131  NA  NA -1.160552775 -0.309205221 -1.548462547 -0.82005932  NA  NA  NA  NA
#> 38   NA  NA  0.219991886  0.298370167  0.129396773          NA  NA  NA  NA  NA
#> 116  NA  NA  1.456308864  1.659854985  0.699267553          NA  NA  NA  NA  NA
#> 189  NA  NA -0.432193146  0.313549164  0.272256971          NA  NA  NA  NA  NA
#> 97   NA  NA  0.579476112 -0.294255255  0.365318802          NA  NA  NA  NA  NA
#> 191  NA  NA  1.230545005  1.882475978           NA          NA  NA  NA  NA  NA
#> 146  NA  NA  2.289545492  0.159415538           NA          NA  NA  NA  NA  NA
#> 165  NA  NA  0.265059463 -0.455800482  1.088920147          NA  NA  NA  NA  NA
#> 185  NA  NA  0.022377309 -0.426522004  0.091857965  0.75956809  NA  NA  NA  NA
#> 138  NA  NA -2.231277272 -0.694795456  0.133839911          NA  NA  NA  NA  NA
#> 153  NA  NA  1.907629641  1.517846450  1.263598570          NA  NA  NA  NA  NA
#> 204  NA  NA  1.061672514  2.650073199  0.151835116          NA  NA  NA  NA  NA
#> 147  NA  NA  0.187860574  0.579368865  0.655907632          NA  NA  NA  NA  NA
#> 125  NA  NA -1.239065746  0.387755522  0.157651335          NA  NA  NA  NA  NA
#> 124  NA  NA -0.456798138 -0.875977109 -1.399862940          NA  NA  NA  NA  NA
#> 199  NA  NA -1.303444865  0.182793562 -0.649631768          NA  NA  NA  NA  NA
#> 170  NA  NA  0.876043121 -0.999617991 -0.731720256          NA  NA  NA  NA  NA
#> 117  NA  NA -1.206198319  0.255426675           NA          NA  NA  NA  NA  NA
#> 175  NA  NA -1.045886695 -1.285392566 -0.775775431          NA  NA  NA  NA  NA
#> 159  NA  NA  0.896952259 -0.615301157 -0.040732817          NA  NA  NA  NA  NA
#> 182  NA  NA -0.586233672  0.944547634           NA          NA  NA  NA  NA  NA
#> 20   NA  NA  1.405337407  1.791978543           NA          NA  NA  NA  NA  NA
#> 181  NA  NA -1.591612789 -0.059487443           NA          NA  NA  NA  NA  NA
#> 109  NA  NA  0.592013647  1.361137445  0.511842761 -0.33116131  NA  NA  NA  NA
#> 78   NA  NA -2.189161083  1.093226115           NA          NA  NA  NA  NA  NA
#> 132  NA  NA -1.223267488 -0.228380799 -1.210741617          NA  NA  NA  NA  NA
#> 56   NA  NA  1.067217042  1.769194054  2.597633378          NA  NA  NA  NA  NA
#> 74   NA  NA  0.707391410  1.822787940           NA          NA  NA  NA  NA  NA
#> 100  NA  NA -0.667157681 -0.040689730 -0.417082933          NA  NA  NA  NA  NA
#> 29   NA  NA  1.654192026  1.203481260           NA          NA  NA  NA  NA  NA
#> 150  NA  NA  0.722132438  0.384198663  1.258952871          NA  NA  NA  NA  NA
#> 203  NA  NA  1.211151785  1.748788434  0.632232829          NA  NA  NA  NA  NA
#> 16   NA  NA -0.392233222  0.580958143 -0.334875341 -1.77542192  NA  NA  NA  NA
#> 183  NA  NA  0.160837124 -0.485648767           NA          NA  NA  NA  NA  NA
#> 85   NA  NA  0.438193135  1.195513529  0.187605703 -0.68127328  NA  NA  NA  NA
#> 135  NA  NA  0.217051162  0.251547284 -0.288067683          NA  NA  NA  NA  NA
#> 171  NA  NA  0.268659115  1.297518339           NA          NA  NA  NA  NA  NA
#> 123  NA  NA  1.650268706  1.108567698           NA          NA  NA  NA  NA  NA
#> 32   NA  NA  0.073408302 -0.051754114           NA          NA  NA  NA  NA  NA
#> 140  NA  NA -0.068888769  0.065824913           NA          NA  NA  NA  NA  NA
#> 112  NA  NA -1.468347821 -1.097619587 -1.131833509          NA  NA  NA  NA  NA
#> 148  NA  NA -0.679230742 -0.769913452           NA          NA  NA  NA  NA  NA
#> 198  NA  NA -0.214727355 -2.068994352 -0.699089584          NA  NA  NA  NA  NA
#> 180  NA  NA  1.582748531  0.678406672           NA          NA  NA  NA  NA  NA
#> 58   NA  NA -0.703449768 -1.730260272 -1.023356087 -1.84073099  NA  NA  NA  NA
#> 177  NA  NA -0.178112964  0.456134599  1.036076815  0.45508508  NA  NA  NA  NA
#>     L_8 L_9        X_2BL        W_2BL        L_2BL  event_time ev   prog_time
#> 73   NA  NA -1.327331255 -0.095487568 -0.794908725 0.002353204  1 0.002353204
#> 89   NA  NA -0.265020949 -1.483969090  1.648630293 0.064847842  1 0.064847842
#> 71   NA  NA -1.129568395 -0.984978713  0.009935887 0.092757381  1 0.092757381
#> 120  NA  NA  0.792347143 -1.736128805  1.251584939 0.196058778  1 0.196058778
#> 77   NA  NA -2.128069870 -2.639347932  0.010827752 0.530073390  1 0.001328669
#> 54   NA  NA  0.485275870 -0.866285096 -1.196253657 0.212843507  1 0.183000212
#> 121  NA  NA -1.177596300  1.333692983 -1.341086424 0.456500274  1 0.456500274
#> 141  NA  NA -0.762932289  0.277604643 -0.556643307 0.734612244  1 0.734612244
#> 70   NA  NA  0.048068709 -1.454541666  1.104915317 0.069884546  1 0.020119427
#> 34   NA  NA -0.036213207 -0.096956372 -0.278771595 0.454292635  0 0.454292635
#> 50   NA  NA -1.624615430  0.749120907 -1.197744279 0.492481324  1 0.015364921
#> 111  NA  NA -1.307237931 -1.390178357  1.058795806 0.188017773  1 0.123881671
#> 94   NA  NA -0.734696342  1.815126136  1.164476460 0.925925952  1 0.745398265
#> 107  NA  NA -0.725479323 -2.203344553  1.480581378 0.890304464  1 0.064051450
#> 19   NA  NA -0.775199253 -0.099402677  0.651809529 0.773393333  0 0.269945919
#> 115  NA  NA -1.230762282  0.245524670  0.677606918 0.635864155  1 0.021325326
#> 23   NA  NA -0.371732264  0.217880795 -0.382772108 1.112603037  1 0.788840848
#> 47   NA  NA -1.543695084 -1.075999051 -0.314709642 0.745791741  1 0.097421106
#> 67   NA  NA -0.417580500  0.002995459 -0.313104798 0.488913226  1 0.488913226
#> 81   NA  NA -0.115655250 -1.231787034  2.094910336 0.413635635  1 0.094178199
#> 142  NA  NA -0.464005889  0.212185450  0.956033786 0.424277238  1 0.274356005
#> 37   NA  NA -0.347298344  0.930313277 -1.677444119 0.680339926  0 0.680339926
#> 8    NA  NA -1.518914931 -0.397932510  1.802232725 0.413270074  1 0.096781811
#> 79   NA  NA -2.299280635  0.181176517 -0.278958962 1.155363670  1 0.055107058
#> 95   NA  NA  0.814690644 -0.681750586 -0.514629206 1.230664890  1 1.230664890
#> 26   NA  NA -0.595717879  0.460046529 -0.969804695 0.800592139  0 0.005766623
#> 155  NA  NA -2.149099081  1.518396696  0.637275130 1.246141018  1 0.589996834
#> 48   NA  NA -0.811621646  0.712021433 -2.270163564 0.396951039  1 0.014544475
#> 158  NA  NA -0.961820003 -2.672681664 -1.477314055 0.285045990  1 0.066308295
#> 200  NA  NA -1.096651989 -0.623431035 -1.155078015 0.873109906  1 0.772258627
#> 134  NA  NA -0.845100201  0.278552309 -0.760119705 0.838567971  1 0.447209453
#> 156  NA  NA -0.169367187  0.382999960 -1.193336784 0.294603671  1 0.294603671
#> 119  NA  NA -0.693481350  0.310601303 -0.043105375 0.240644112  0 0.240644112
#> 127  NA  NA -1.775397920  0.331400128  1.441918671 1.499342072  1 0.190921027
#> 103  NA  NA  0.415851412  0.163243552  0.660727284 0.097913614  0 0.097913614
#> 193  NA  NA  0.210213656 -0.700733592 -0.874565020 1.050205320  0 0.065140308
#> 118  NA  NA  0.132161461 -0.399372778 -1.821175364 0.392917264  1 0.392917264
#> 209  NA  NA -1.534128183 -0.220433050  0.151799854 0.774877648  1 0.774877648
#> 98   NA  NA -0.520123482 -0.363756203 -0.533248173 0.366329872  1 0.081656717
#> 61   NA  NA -1.793427835 -1.149284301  2.128569801 0.312165626  1 0.088841900
#> 4    NA  NA  1.920388895 -1.435639510  0.818440069 0.322366037  0 0.322366037
#> 157  NA  NA -0.829284774 -2.197902005  0.177990132 0.401170359  1 0.241695793
#> 7    NA  NA -0.803455568 -1.156774088  1.141601107 0.632711979  1 0.632711979
#> 66   NA  NA -2.696984957 -0.649366033  0.098571479 0.524864219  1 0.040851099
#> 60   NA  NA -0.657858374  1.193529373 -1.132853826 1.167978394  1 0.550992077
#> 43   NA  NA  0.226715318 -2.325448321 -0.053452190 1.202545810  1 0.028023311
#> 210  NA  NA  1.936517565 -1.682346208  0.438172274 0.785936024  1 0.785936024
#> 28   NA  NA  0.677246175 -0.390359715  0.854217903 1.199394761  1 1.199394761
#> 17   NA  NA -0.711604050  0.114119709  0.222469688 1.956451685  1 0.852180502
#> 162  NA  NA -0.303133823 -1.762610526  0.645166675 0.483277335  1 0.356613289
#> 205  NA  NA  0.743963173  0.165505083  0.667571054 0.984376274  1 0.984376274
#> 80   NA  NA -0.631220028 -0.659518318 -1.498630983 0.311257623  1 0.311257623
#> 188  NA  NA  0.860846390  0.112258388  1.274302001 1.267754245  1 1.267754245
#> 25   NA  NA  0.039356836 -0.729974507  0.020275989 0.423589579  1 0.239496524
#> 84   NA  NA -1.291981751  0.512433788 -0.659231022 1.607915078  1 0.028150883
#> 108  NA  NA -0.523347405  0.256595324  0.571464822 2.004948019  1 0.703012384
#> 9    NA  NA -1.555776031 -1.606157999 -1.759947009 0.550429530  1 0.044363078
#> 144  NA  NA -0.213436234 -2.153929733  0.902280636 0.503631211  1 0.503631211
#> 195  NA  NA  0.830192409 -0.066677906  0.495780891 1.353885638  0 0.787881073
#> 168  NA  NA  0.126316649  0.289605229 -0.763095428 1.991738812  1 0.348125722
#> 149  NA  NA -0.671788473  0.228321508  0.111138463 2.558939626  1 0.026259135
#> 41   NA  NA -0.775025328  1.029129162  0.512247411 1.636406122  0 0.207894537
#> 113  NA  NA  0.164531916  0.534821963 -0.325855498 1.533094289  1 0.781922153
#> 167  NA  NA -0.325262792  0.358871420  0.689875751 2.362955308  1 0.280843441
#> 44   NA  NA -0.711380049  0.667368499 -1.126027190 1.966324345  0 0.030605681
#> 166  NA  NA -1.425972572 -1.665574147 -0.464071382 1.558084983  1 0.040692388
#> 55   NA  NA -0.132628455 -0.472669441  0.153955408 1.158489489  1 0.911210617
#> 92   NA  NA  0.532975205 -1.462279621  0.352098754 1.111758823  1 1.111758823
#> 197  NA  NA  0.698499194  0.678753415  0.478430279 2.654796504  0 2.654796504
#> 143  NA  NA -1.950884497  2.074215381  1.297360843 2.111744753  0 0.205602082
#> 12   NA  NA -1.343792501 -0.017609808  0.180525062 1.915677032  1 0.539439651
#> 164  NA  NA  0.301795118  0.151650132 -0.753162179 2.973012587  0 2.221551469
#> 22   NA  NA -0.299127703 -0.955102753 -0.145914654 1.477337829  1 0.215601515
#> 42   NA  NA -0.153702393 -1.553093354  1.090507677 2.001340019  1 1.253366427
#> 46   NA  NA  0.084716934  0.939843814  1.456810877 1.458892201  1 0.396798704
#> 68   NA  NA  0.175446188 -0.722859130  1.318761959 1.691990624  1 0.833632933
#> 18   NA  NA  1.005909127 -0.069511408  1.224726357 1.427555383  1 1.400994310
#> 179  NA  NA  0.104606714  1.151982609  0.927508093 2.914007666  1 2.642649693
#> 75   NA  NA  1.086098823  0.365714231 -0.357133375 2.032890409  1 2.032890409
#> 96   NA  NA  0.629507845 -0.606809943 -0.562968455 3.039800692  1 0.447685385
#> 49   NA  NA -0.223436252  0.039857822 -0.072309982 1.792127442  0 0.635828158
#> 33   NA  NA -0.232978661  0.214883569  0.116506419 2.095950320  0 0.827211563
#> 129  NA  NA -1.208932144  1.000724557  0.123690580 1.692466029  0 0.307234605
#> 122  NA  NA  0.194159505 -2.603336740  1.014189038 1.965424873  0 0.662039293
#> 102  NA  NA -0.704295442  0.854427448  0.402181441 3.298269080  0 2.051598582
#> 51   NA  NA -1.298458855  0.272490237 -1.822887636 3.032966860  0 1.117471991
#> 88   NA  NA  0.823198729 -1.425318914 -1.317752176 1.940475225  0 0.252115286
#> 161  NA  NA  0.521056693  0.493914632  2.320931135 1.577465885  0 1.141087445
#> 160  NA  NA -0.784471147 -0.122231761 -0.601188181 2.938799294  0 2.593368590
#> 105  NA  NA  0.898586893 -0.792793107 -0.272388414 1.559098005  0 0.697569114
#> 86   NA  NA -1.189930229 -0.247629059 -0.222495148 1.587660448  0 0.072223624
#> 10   NA  NA -0.534201600  0.469382921 -0.147538958 3.036745627  0 0.633187685
#> 15   NA  NA -1.790069456  0.527380537 -0.017123892 2.404734716  0 0.155656004
#> 93   NA  NA -0.202056466 -2.093849242 -0.843192288 3.146626260  0 3.146626260
#> 64   NA  NA  1.308854820  0.327536272 -0.475192101 2.065269654  0 1.503553151
#> 110  NA  NA -1.744183492 -0.009079277 -1.652279469 3.289224584  0 0.391421074
#> 62   NA  NA -1.033135468  1.026890967 -0.633371120 1.976171788  0 0.195317891
#> 172  NA  NA  0.290301415  0.570356577  0.494473021 3.392820178  0 2.151593955
#> 174  NA  NA -0.834203453  0.847661951 -1.302929075 2.416213935  0 0.643556867
#> 192  NA  NA  0.572901813  0.763149796 -1.526714454 2.691348348  0 0.398635840
#> 106  NA  NA  1.332554213  0.478843397 -0.592491787 2.699221319  0 0.864803851
#> 82   NA  NA  0.690667743  0.065509091 -1.266096352 3.236655228  0 0.240744960
#> 190  NA  NA -1.745148998 -0.842275576 -0.827964487 3.131688696  0 0.418314372
#> 39   NA  NA -1.151406159  1.617593486  0.446124915 3.199698205  0 0.031551863
#> 45   NA  NA -0.473650038 -0.933266617  0.066932290 3.072801209  0 0.015857148
#> 87   NA  NA -0.364963550  0.080245600  0.761639720 2.739771616  0 0.012504875
#> 30   NA  NA -1.150449433  0.291347371  0.944487137 2.253803359  0 0.221457293
#> 69   NA  NA  0.447675407 -1.903557577 -1.956659025 2.109129844  0 0.103180898
#> 169  NA  NA -0.066396775 -0.904853451 -0.519461812 3.070887534  0 0.378734595
#> 202  NA  NA  0.574880855 -0.350173032  2.028201820 2.615857600  0 0.750348555
#> 196  NA  NA -0.391877880  0.926434516  0.980772093 1.832194052  0 0.973348224
#> 133  NA  NA -0.235937890 -0.373016726  1.381269804 2.710601567  0 1.706690859
#> 2    NA  NA -0.618695501 -0.772981776  0.646887889 3.439330093  0 1.256963971
#> 186  NA  NA  0.607833809  0.705479091 -0.862364263 2.316735843  0 0.886289405
#> 208  NA  NA -1.930423709  0.828217986  0.228503882 2.194600088  0 0.033820246
#> 35   NA  NA -0.172983511 -1.146577802 -0.760092491 1.487897399  0 0.681201926
#> 163  NA  NA  0.250827790  0.002250473 -0.682565693 1.730919743  0 1.730919743
#> 31   NA  NA -1.162498817 -1.041109234  1.064098669 1.491037920  0 0.036723780
#> 114  NA  NA -0.690935395 -0.853781551  0.452656967 1.988938922  0 0.551298673
#> 57   NA  NA -0.234955033  0.513090275 -1.625695619 3.009433361  0 2.705799648
#> 3    NA  NA -0.657754469 -0.818070273  0.782547913 2.282342987  0 0.180068310
#> 184  NA  NA -0.216609744 -1.179172756 -0.736173249 2.526820653  0 0.859066252
#> 24   NA  NA -1.987654126  0.870014181  0.501218212 1.997173630  0 0.204938749
#> 76   NA  NA -0.354982605  0.405993998  2.060636237 3.108176357  0 0.091543447
#> 63   NA  NA  0.848124435  0.122788989 -1.029007487 2.593301453  0 1.210957521
#> 11   NA  NA  0.868720605 -1.340219960  0.917805322 2.254999778  0 0.640899669
#> 101  NA  NA -0.132253303 -0.960536609 -0.788186265 2.448123648  0 1.742991241
#> 13   NA  NA -1.107052754 -0.160108356 -0.531805789 2.375588767  0 1.009914703
#> 201  NA  NA  0.472093208 -2.315046633 -2.182697279 1.624113011  0 0.051545083
#> 152  NA  NA -0.228282333  0.657808273  0.139604637 2.745310251  0 2.011334908
#> 27   NA  NA -1.606942139  1.620908680 -0.084679412 2.930390294  0 0.674797021
#> 1    NA  NA -0.223839052  1.451126910 -0.949173836 2.639211520  0 2.639211520
#> 6    NA  NA  1.239778480  1.469164347  0.084687426 3.213883845  0 0.275590032
#> 59   NA  NA  1.190271478  1.903628559 -1.015170822 3.225515027  0 3.062311968
#> 90   NA  NA  1.375458080 -1.046695306  0.374152164 2.813407329  0 1.081900943
#> 176  NA  NA  0.086903103 -0.353964586 -2.820126313 3.117524430  0 0.206618948
#> 40   NA  NA  0.367930480 -0.186063990  0.319532233 1.662267254  0 1.214121793
#> 99   NA  NA  0.182612728 -0.540114035 -0.046764494 1.889753358  0 0.185104044
#> 14   NA  NA  0.391558482 -0.886896115  1.397524883 2.195999028  0 0.133196334
#> 126  NA  NA  0.748941332 -0.047750442  1.835764370 1.992640784  0 0.372663921
#> 206  NA  NA -0.263404679 -0.866087336  1.662589287 2.595806142  0 0.395334046
#> 104  NA  NA  0.402476845 -1.013155029  0.745789066 2.648999685  0 0.072042659
#> 154  NA  NA -0.769931534  0.574902490 -0.740002546 1.814299107  0 0.051951131
#> 151  NA  NA -0.430827986 -0.038277348 -1.242959721 2.194467321  0 2.163924950
#> 139  NA  NA  0.965517293 -0.784660120 -0.651631175 2.908149838  0 2.908149838
#> 91   NA  NA  2.936759868 -1.854374754 -1.196621759 2.754422292  0 0.056830628
#> 36   NA  NA -0.539366090 -0.941360999  1.915441266 1.502127951  0 0.037793023
#> 21   NA  NA -0.093649019 -0.026483949 -0.388766460 1.815228232  0 0.998794142
#> 136  NA  NA  1.291999643  0.726023557  0.080080154 1.828098366  0 1.664909047
#> 53   NA  NA  0.483472595  0.228799037 -0.378212300 1.533185803  0 0.844724129
#> 187  NA  NA  0.000239312  0.108797712 -0.904551243 1.674409366  0 1.674409366
#> 207  NA  NA -0.419221617 -0.201076025  0.492436120 3.450517944  0 1.648361543
#> 52   NA  NA  0.362321932  0.836470166  0.339400827 2.609274742  0 1.596512601
#> 173  NA  NA  0.329077694 -1.424230301 -0.323890827 2.933088858  0 0.028917838
#> 83   NA  NA  1.539132115 -1.317665756  0.759673191 3.143345806  0 0.104111469
#> 72   NA  NA -0.025271504 -0.890257399  0.945645245 1.867180959  0 0.257842064
#> 211  NA  NA -0.889314955  1.085570464 -0.359191132 2.003891013  0 0.608609676
#> 65   NA  NA -1.058719062  0.086101145 -2.120953225 2.933514296  0 0.180124503
#> 130  NA  NA  0.002883867  1.050134106 -1.531253617 2.382442025  0 1.012508577
#> 128  NA  NA -0.338679322  0.742394690  2.993528520 1.700275874  0 1.296496445
#> 178  NA  NA -0.013440081  0.520021399  0.218882485 2.886291088  0 0.380344274
#> 145  NA  NA -0.605216897  0.216480811 -0.784392191 1.966119197  0 0.591917714
#> 137  NA  NA  0.157472272  0.429029574  1.424297179 2.339596255  0 1.884994603
#> 194  NA  NA -0.042452885 -0.642092390 -1.531818574 1.519077732  0 1.519077732
#> 5    NA  NA -1.188011877  0.507616069  0.892913949 2.638203583  0 0.651173859
#> 131  NA  NA  1.283452883  1.329820908  0.599701684 3.356665129  0 3.356665129
#> 38   NA  NA  1.024158477  1.407212112 -0.453391305 2.094716181  0 2.094716181
#> 116  NA  NA  0.532525947 -0.620634085  0.844516876 2.912181630  0 2.912181630
#> 189  NA  NA -0.651318395  1.910394287 -0.432193146 2.813197410  0 0.068286222
#> 97   NA  NA  2.732720118 -1.764973150  0.365318802 2.552347922  0 2.552347922
#> 191  NA  NA  0.190260773  0.292436501  1.678640978 1.952491022  0 1.952491022
#> 146  NA  NA  0.495023102  1.497935962  2.289545492 1.888812617  0 0.957451626
#> 165  NA  NA  0.291779038  0.548449812  0.265059463 2.879067909  0 0.095055484
#> 185  NA  NA -0.278773118  1.086827859 -0.426522004 3.374368395  0 1.432218767
#> 138  NA  NA -0.555270884  1.327477435 -0.694795456 2.548222541  0 1.145872570
#> 153  NA  NA  0.422243636  1.002425646  1.907629641 2.020404576  0 0.566880720
#> 204  NA  NA -1.138775086  3.122498290  0.460231707 2.355819957  0 2.355819957
#> 147  NA  NA  0.051521241 -0.348442242  0.187860574 2.063621580  0 0.691861665
#> 125  NA  NA -0.334662851  0.374995053  0.157651335 2.130415711  0 2.130415711
#> 124  NA  NA  0.637578012  1.689430313 -2.377902814 2.364459189  0 2.364459189
#> 199  NA  NA  0.219953958 -2.303321618 -1.303444865 2.915514640  0 0.868414115
#> 170  NA  NA  1.070618491 -0.824807481 -1.725631310 2.018080134  0 2.018080134
#> 117  NA  NA -0.041284567  0.108344013 -1.206198319 1.643331278  0 0.047817689
#> 175  NA  NA  0.194415453  0.054490877 -1.045886695 2.417483228  0 0.685151471
#> 159  NA  NA  1.982080939 -1.564794045  0.743934168 2.492320854  0 2.492320854
#> 182  NA  NA -0.802574794  0.516508078 -0.586233672 1.754409486  0 0.952706289
#> 20   NA  NA  0.677665541 -0.137006008  1.791978543 1.848034045  0 1.828974557
#> 181  NA  NA  0.944153820  0.075108124 -1.591612789 1.610859601  0 0.612705867
#> 109  NA  NA -0.461941033  0.412373712  0.592013647 3.048073172  0 0.689042420
#> 78   NA  NA  0.768800290 -1.727281923  0.197220746 1.949730885  0 1.949730885
#> 132  NA  NA  0.193206155 -0.573606808 -1.223267488 2.123299163  0 0.489980717
#> 56   NA  NA  0.306189764 -0.021299324  1.067217042 2.693983021  0 0.891280987
#> 74   NA  NA  2.224260498 -0.535121224  0.778350480 1.937423018  0 1.937423018
#> 100  NA  NA  0.262549454 -0.534797772 -0.667157681 2.234270757  0 0.838104308
#> 29   NA  NA  0.594149748 -1.154317933  1.654192026 1.734116142  0 0.094964164
#> 150  NA  NA  1.065468793  0.828604184  0.384198663 2.624894429  0 1.896446076
#> 203  NA  NA  1.707458672 -0.222328861  0.306429525 2.870742497  0 2.870742497
#> 16   NA  NA  0.294382832 -0.118895848 -0.392233222 3.463437349  0 0.710660583
#> 183  NA  NA -0.529701796 -0.425578710  0.160837124 1.547411234  0 0.342474712
#> 85   NA  NA -0.236723170 -0.481200350  0.438193135 3.043860347  0 0.419877629
#> 135  NA  NA  1.034367385 -1.229962916  0.217051162 2.714988536  0 0.360640278
#> 171  NA  NA  1.043612025 -0.555775125  1.297518339 1.556319889  0 1.556319889
#> 123  NA  NA  1.330582986  0.358559697  2.532354356 1.853880733  0 1.853880733
#> 32   NA  NA  0.432956814  0.564777186 -0.844519335 1.983068798  0 1.983068798
#> 140  NA  NA -1.229867169 -0.472478384 -0.068888769 1.665965411  0 0.180869230
#> 112  NA  NA  0.133645058 -1.249604658 -1.468347821 2.242371005  0 0.090253978
#> 148  NA  NA  3.108110657  0.655518907 -1.028741210 1.685644908  0 1.685644908
#> 198  NA  NA  0.058599917  0.102162280 -0.214727355 2.701082746  0 0.771312046
#> 180  NA  NA  2.954902872  0.876522518  1.011664201 1.881428289  0 1.881428289
#> 58   NA  NA  1.630736201  0.860095869 -0.703449768 3.258573407  0 0.636626512
#> 177  NA  NA  0.817882922  0.495241738  0.552667165 3.203632067  0 3.203632067
#>     prog_ev switch admin_cens calendar_start_time calendar_end_of_study
#> 73        0      0      FALSE         0.175900844              3.469317
#> 89        0      0      FALSE         0.240610989              3.469317
#> 71        0      0      FALSE         0.310362298              3.469317
#> 120       0      0      FALSE         0.230675962              3.469317
#> 77        1      0      FALSE         0.129561043              3.469317
#> 54        1      1      FALSE         0.486544588              3.469317
#> 121       0      0      FALSE         0.258812577              3.469317
#> 141       0      0      FALSE         0.017892577              3.469317
#> 70        1      0      FALSE         0.781446745              3.469317
#> 34        0      0      FALSE         0.409292586              3.469317
#> 50        1      0      FALSE         0.379280732              3.469317
#> 111       1      0      FALSE         0.696683944              3.469317
#> 94        1      1      FALSE         0.119390312              3.469317
#> 107       1      0      FALSE         0.167973878              3.469317
#> 19        1      1      FALSE         0.287791566              3.469317
#> 115       1      0      FALSE         0.456148904              3.469317
#> 23        1      0      FALSE         0.027233258              3.469317
#> 47        1      0      FALSE         0.398542133              3.469317
#> 67        0      0      FALSE         0.690654131              3.469317
#> 81        1      0      FALSE         0.788881110              3.469317
#> 142       1      0      FALSE         0.788018044              3.469317
#> 37        0      0      FALSE         0.555874887              3.469317
#> 8         1      0      FALSE         0.912976371              3.469317
#> 79        1      0      FALSE         0.209692389              3.469317
#> 95        0      0      FALSE         0.153109332              3.469317
#> 26        1      1      FALSE         0.634865209              3.469317
#> 155       1      0      FALSE         0.201184371              3.469317
#> 48        1      0      FALSE         1.117846210              3.469317
#> 158       1      0      FALSE         1.259704937              3.469317
#> 200       1      0      FALSE         0.673012327              3.469317
#> 134       1      0      FALSE         0.707699491              3.469317
#> 156       0      0      FALSE         1.425408739              3.469317
#> 119       0      0      FALSE         1.526917506              3.469317
#> 127       1      0      FALSE         0.290689775              3.469317
#> 103       0      0      FALSE         1.739278159              3.469317
#> 193       1      0      FALSE         0.794846342              3.469317
#> 118       0      0      FALSE         1.476313948              3.469317
#> 209       0      0      FALSE         1.108871426              3.469317
#> 98        1      0      FALSE         1.545478593              3.469317
#> 61        1      0      FALSE         1.620034405              3.469317
#> 4         0      0      FALSE         1.683585772              3.469317
#> 157       1      0      FALSE         1.630332469              3.469317
#> 7         0      0      FALSE         1.460077809              3.469317
#> 66        1      0      FALSE         1.570396788              3.469317
#> 60        1      0      FALSE         0.967974691              3.469317
#> 43        1      0      FALSE         0.954967133              3.469317
#> 210       0      0      FALSE         1.376560731              3.469317
#> 28        0      0      FALSE         0.973547168              3.469317
#> 17        1      1      FALSE         0.224437499              3.469317
#> 162       1      0      FALSE         1.763638028              3.469317
#> 205       0      0      FALSE         1.263558946              3.469317
#> 80        0      0      FALSE         1.962192562              3.469317
#> 188       0      0      FALSE         1.011096295              3.469317
#> 25        1      1      FALSE         1.931264218              3.469317
#> 84        1      1      FALSE         0.766086280              3.469317
#> 108       1      0      FALSE         0.391459564              3.469317
#> 9         1      0      FALSE         1.859801196              3.469317
#> 144       0      0      FALSE         1.966002672              3.469317
#> 195       1      0      FALSE         1.196368776              3.469317
#> 168       1      0      FALSE         0.648998660              3.469317
#> 149       1      0      FALSE         0.097755728              3.469317
#> 41        1      0      FALSE         1.021000893              3.469317
#> 113       1      0      FALSE         1.187621335              3.469317
#> 167       1      0      FALSE         0.471971720              3.469317
#> 44        1      0      FALSE         0.888675342              3.469317
#> 166       1      0      FALSE         1.298563391              3.469317
#> 55        1      1      FALSE         1.699153261              3.469317
#> 92        0      0      FALSE         1.852669914              3.469317
#> 197       0      0      FALSE         0.330281170              3.469317
#> 143       1      0      FALSE         0.882662911              3.469317
#> 12        1      0      FALSE         1.153891793              3.469317
#> 164       1      0      FALSE         0.099831816              3.469317
#> 22        1      0      FALSE         1.639164208              3.469317
#> 42        1      1      FALSE         1.165133607              3.469317
#> 46        1      0      FALSE         1.721893364              3.469317
#> 68        1      0      FALSE         1.491610168              3.469317
#> 18        1      0      FALSE         1.773076950              3.469317
#> 179       1      0      FALSE         0.422802716              3.469317
#> 75        0      0      FALSE         1.339436167              3.469317
#> 96        1      0      FALSE         0.429516181              3.469317
#> 49        1      1       TRUE         1.677189431              3.469317
#> 33        1      1       TRUE         1.373366553              3.469317
#> 129       1      0       TRUE         1.776850844              3.469317
#> 122       1      0       TRUE         1.503892000              3.469317
#> 102       1      0       TRUE         0.171047793              3.469317
#> 51        1      1       TRUE         0.436350013              3.469317
#> 88        1      0       TRUE         1.528841648              3.469317
#> 161       1      0       TRUE         1.891850988              3.469317
#> 160       1      0       TRUE         0.530517579              3.469317
#> 105       1      0       TRUE         1.910218868              3.469317
#> 86        1      0       TRUE         1.881656425              3.469317
#> 10        1      1       TRUE         0.432571246              3.469317
#> 15        1      0       TRUE         1.064582157              3.469317
#> 93        0      0       TRUE         0.322690613              3.469317
#> 64        1      0       TRUE         1.404047219              3.469317
#> 110       1      0       TRUE         0.180092289              3.469317
#> 62        1      0       TRUE         1.493145085              3.469317
#> 172       1      0       TRUE         0.076496695              3.469317
#> 174       1      0       TRUE         1.053102938              3.469317
#> 192       1      0       TRUE         0.777968525              3.469317
#> 106       1      0       TRUE         0.770095554              3.469317
#> 82        1      1       TRUE         0.232661645              3.469317
#> 190       1      0       TRUE         0.337628177              3.469317
#> 39        1      1       TRUE         0.269618668              3.469317
#> 45        1      1       TRUE         0.396515664              3.469317
#> 87        1      1       TRUE         0.729545257              3.469317
#> 30        1      0       TRUE         1.215513514              3.469317
#> 69        1      0       TRUE         1.360187029              3.469317
#> 169       1      0       TRUE         0.398429339              3.469317
#> 202       1      0       TRUE         0.853459273              3.469317
#> 196       1      0       TRUE         1.637122821              3.469317
#> 133       1      0       TRUE         0.758715306              3.469317
#> 2         1      1       TRUE         0.029986780              3.469317
#> 186       1      0       TRUE         1.152581030              3.469317
#> 208       1      0       TRUE         1.274716785              3.469317
#> 35        1      0       TRUE         1.981419474              3.469317
#> 163       0      0       TRUE         1.738397130              3.469317
#> 31        1      0       TRUE         1.978278953              3.469317
#> 114       1      0       TRUE         1.480377951              3.469317
#> 57        1      1       TRUE         0.459883512              3.469317
#> 3         1      1       TRUE         1.186973886              3.469317
#> 184       1      0       TRUE         0.942496220              3.469317
#> 24        1      1       TRUE         1.472143243              3.469317
#> 76        1      1       TRUE         0.361140516              3.469317
#> 63        1      0       TRUE         0.876015420              3.469317
#> 11        1      0       TRUE         1.214317095              3.469317
#> 101       1      0       TRUE         1.021193225              3.469317
#> 13        1      1       TRUE         1.093728106              3.469317
#> 201       1      0       TRUE         1.845203862              3.469317
#> 152       1      0       TRUE         0.724006622              3.469317
#> 27        1      1       TRUE         0.538926579              3.469317
#> 1         0      0       TRUE         0.830105353              3.469317
#> 6         1      1       TRUE         0.255433028              3.469317
#> 59        1      1       TRUE         0.243801846              3.469317
#> 90        1      0       TRUE         0.655909544              3.469317
#> 176       1      0       TRUE         0.351792443              3.469317
#> 40        1      0       TRUE         1.807049619              3.469317
#> 99        1      0       TRUE         1.579563515              3.469317
#> 14        1      0       TRUE         1.273317845              3.469317
#> 126       1      0       TRUE         1.476676089              3.469317
#> 206       1      0       TRUE         0.873510731              3.469317
#> 104       1      0       TRUE         0.820317188              3.469317
#> 154       1      0       TRUE         1.655017766              3.469317
#> 151       1      0       TRUE         1.274849552              3.469317
#> 139       0      0       TRUE         0.561167035              3.469317
#> 91        1      0       TRUE         0.714894581              3.469317
#> 36        1      1       TRUE         1.967188922              3.469317
#> 21        1      1       TRUE         1.654088641              3.469317
#> 136       1      0       TRUE         1.641218507              3.469317
#> 53        1      1       TRUE         1.936131070              3.469317
#> 187       0      0       TRUE         1.794907507              3.469317
#> 207       1      0       TRUE         0.018798929              3.469317
#> 52        1      1       TRUE         0.860042131              3.469317
#> 173       1      0       TRUE         0.536228015              3.469317
#> 83        1      1       TRUE         0.325971067              3.469317
#> 72        1      0       TRUE         1.602135914              3.469317
#> 211       1      0       TRUE         1.465425860              3.469317
#> 65        1      1       TRUE         0.535802577              3.469317
#> 130       1      0       TRUE         1.086874848              3.469317
#> 128       1      0       TRUE         1.769040999              3.469317
#> 178       1      0       TRUE         0.583025785              3.469317
#> 145       1      0       TRUE         1.503197676              3.469317
#> 137       1      0       TRUE         1.129720618              3.469317
#> 194       0      0       TRUE         1.950239141              3.469317
#> 5         1      1       TRUE         0.831113290              3.469317
#> 131       0      0       TRUE         0.112651744              3.469317
#> 38        0      0       TRUE         1.374600692              3.469317
#> 116       0      0       TRUE         0.557135243              3.469317
#> 189       1      0       TRUE         0.656119463              3.469317
#> 97        0      0       TRUE         0.916968951              3.469317
#> 191       0      0       TRUE         1.516825851              3.469317
#> 146       1      0       TRUE         1.580504256              3.469317
#> 165       1      0       TRUE         0.590248964              3.469317
#> 185       1      0       TRUE         0.094948479              3.469317
#> 138       1      0       TRUE         0.921094332              3.469317
#> 153       1      0       TRUE         1.448912297              3.469317
#> 204       0      0       TRUE         1.113496916              3.469317
#> 147       1      0       TRUE         1.405695293              3.469317
#> 125       0      0       TRUE         1.338901162              3.469317
#> 124       0      0       TRUE         1.104857684              3.469317
#> 199       1      0       TRUE         0.553802233              3.469317
#> 170       0      0       TRUE         1.451236739              3.469317
#> 117       1      0       TRUE         1.825985595              3.469317
#> 175       1      0       TRUE         1.051833645              3.469317
#> 159       0      0       TRUE         0.976996019              3.469317
#> 182       1      0       TRUE         1.714907387              3.469317
#> 20        1      0       TRUE         1.621282828              3.469317
#> 181       1      0       TRUE         1.858457272              3.469317
#> 109       1      0       TRUE         0.421243701              3.469317
#> 78        0      0       TRUE         1.519585988              3.469317
#> 132       1      0       TRUE         1.346017710              3.469317
#> 56        1      1       TRUE         0.775333852              3.469317
#> 74        0      0       TRUE         1.531893855              3.469317
#> 100       1      0       TRUE         1.235046116              3.469317
#> 29        1      1       TRUE         1.735200731              3.469317
#> 150       1      0       TRUE         0.844422444              3.469317
#> 203       0      0       TRUE         0.598574376              3.469317
#> 16        1      1       TRUE         0.005879524              3.469317
#> 183       1      0       TRUE         1.921905639              3.469317
#> 85        1      0       TRUE         0.425456526              3.469317
#> 135       1      0       TRUE         0.754328337              3.469317
#> 171       0      0       TRUE         1.912996984              3.469317
#> 123       0      0       TRUE         1.615436140              3.469317
#> 32        0      0       TRUE         1.486248075              3.469317
#> 140       1      0       TRUE         1.803351462              3.469317
#> 112       1      0       TRUE         1.226945868              3.469317
#> 148       0      0       TRUE         1.783671965              3.469317
#> 198       1      0       TRUE         0.768234127              3.469317
#> 180       0      0       TRUE         1.587888584              3.469317
#> 58        1      1       TRUE         0.210743466              3.469317
#> 177       0      0       TRUE         0.265684806              3.469317
#>     random_cens ev_soll ev_obs
#> 73        FALSE      66     66
#> 89        FALSE      66     66
#> 71        FALSE      66     66
#> 120       FALSE      66     66
#> 77        FALSE      66     66
#> 54        FALSE      66     66
#> 121       FALSE      66     66
#> 141       FALSE      66     66
#> 70        FALSE      66     66
#> 34         TRUE      66     66
#> 50        FALSE      66     66
#> 111       FALSE      66     66
#> 94        FALSE      66     66
#> 107       FALSE      66     66
#> 19         TRUE      66     66
#> 115       FALSE      66     66
#> 23        FALSE      66     66
#> 47        FALSE      66     66
#> 67        FALSE      66     66
#> 81        FALSE      66     66
#> 142       FALSE      66     66
#> 37         TRUE      66     66
#> 8         FALSE      66     66
#> 79        FALSE      66     66
#> 95        FALSE      66     66
#> 26         TRUE      66     66
#> 155       FALSE      66     66
#> 48        FALSE      66     66
#> 158       FALSE      66     66
#> 200       FALSE      66     66
#> 134       FALSE      66     66
#> 156       FALSE      66     66
#> 119        TRUE      66     66
#> 127       FALSE      66     66
#> 103        TRUE      66     66
#> 193        TRUE      66     66
#> 118       FALSE      66     66
#> 209       FALSE      66     66
#> 98        FALSE      66     66
#> 61        FALSE      66     66
#> 4          TRUE      66     66
#> 157       FALSE      66     66
#> 7         FALSE      66     66
#> 66        FALSE      66     66
#> 60        FALSE      66     66
#> 43        FALSE      66     66
#> 210       FALSE      66     66
#> 28        FALSE      66     66
#> 17        FALSE      66     66
#> 162       FALSE      66     66
#> 205       FALSE      66     66
#> 80        FALSE      66     66
#> 188       FALSE      66     66
#> 25        FALSE      66     66
#> 84        FALSE      66     66
#> 108       FALSE      66     66
#> 9         FALSE      66     66
#> 144       FALSE      66     66
#> 195        TRUE      66     66
#> 168       FALSE      66     66
#> 149       FALSE      66     66
#> 41         TRUE      66     66
#> 113       FALSE      66     66
#> 167       FALSE      66     66
#> 44         TRUE      66     66
#> 166       FALSE      66     66
#> 55        FALSE      66     66
#> 92        FALSE      66     66
#> 197        TRUE      66     66
#> 143        TRUE      66     66
#> 12        FALSE      66     66
#> 164        TRUE      66     66
#> 22        FALSE      66     66
#> 42        FALSE      66     66
#> 46        FALSE      66     66
#> 68        FALSE      66     66
#> 18        FALSE      66     66
#> 179       FALSE      66     66
#> 75        FALSE      66     66
#> 96        FALSE      66     66
#> 49        FALSE      66     66
#> 33        FALSE      66     66
#> 129       FALSE      66     66
#> 122       FALSE      66     66
#> 102       FALSE      66     66
#> 51        FALSE      66     66
#> 88        FALSE      66     66
#> 161       FALSE      66     66
#> 160       FALSE      66     66
#> 105       FALSE      66     66
#> 86        FALSE      66     66
#> 10        FALSE      66     66
#> 15        FALSE      66     66
#> 93        FALSE      66     66
#> 64        FALSE      66     66
#> 110       FALSE      66     66
#> 62        FALSE      66     66
#> 172       FALSE      66     66
#> 174       FALSE      66     66
#> 192       FALSE      66     66
#> 106       FALSE      66     66
#> 82        FALSE      66     66
#> 190       FALSE      66     66
#> 39        FALSE      66     66
#> 45        FALSE      66     66
#> 87        FALSE      66     66
#> 30        FALSE      66     66
#> 69        FALSE      66     66
#> 169       FALSE      66     66
#> 202       FALSE      66     66
#> 196       FALSE      66     66
#> 133       FALSE      66     66
#> 2         FALSE      66     66
#> 186       FALSE      66     66
#> 208       FALSE      66     66
#> 35        FALSE      66     66
#> 163       FALSE      66     66
#> 31        FALSE      66     66
#> 114       FALSE      66     66
#> 57        FALSE      66     66
#> 3         FALSE      66     66
#> 184       FALSE      66     66
#> 24        FALSE      66     66
#> 76        FALSE      66     66
#> 63        FALSE      66     66
#> 11        FALSE      66     66
#> 101       FALSE      66     66
#> 13        FALSE      66     66
#> 201       FALSE      66     66
#> 152       FALSE      66     66
#> 27        FALSE      66     66
#> 1         FALSE      66     66
#> 6         FALSE      66     66
#> 59        FALSE      66     66
#> 90        FALSE      66     66
#> 176       FALSE      66     66
#> 40        FALSE      66     66
#> 99        FALSE      66     66
#> 14        FALSE      66     66
#> 126       FALSE      66     66
#> 206       FALSE      66     66
#> 104       FALSE      66     66
#> 154       FALSE      66     66
#> 151       FALSE      66     66
#> 139       FALSE      66     66
#> 91        FALSE      66     66
#> 36        FALSE      66     66
#> 21        FALSE      66     66
#> 136       FALSE      66     66
#> 53        FALSE      66     66
#> 187       FALSE      66     66
#> 207       FALSE      66     66
#> 52        FALSE      66     66
#> 173       FALSE      66     66
#> 83        FALSE      66     66
#> 72        FALSE      66     66
#> 211       FALSE      66     66
#> 65        FALSE      66     66
#> 130       FALSE      66     66
#> 128       FALSE      66     66
#> 178       FALSE      66     66
#> 145       FALSE      66     66
#> 137       FALSE      66     66
#> 194       FALSE      66     66
#> 5         FALSE      66     66
#> 131       FALSE      66     66
#> 38        FALSE      66     66
#> 116       FALSE      66     66
#> 189       FALSE      66     66
#> 97        FALSE      66     66
#> 191       FALSE      66     66
#> 146       FALSE      66     66
#> 165       FALSE      66     66
#> 185       FALSE      66     66
#> 138       FALSE      66     66
#> 153       FALSE      66     66
#> 204       FALSE      66     66
#> 147       FALSE      66     66
#> 125       FALSE      66     66
#> 124       FALSE      66     66
#> 199       FALSE      66     66
#> 170       FALSE      66     66
#> 117       FALSE      66     66
#> 175       FALSE      66     66
#> 159       FALSE      66     66
#> 182       FALSE      66     66
#> 20        FALSE      66     66
#> 181       FALSE      66     66
#> 109       FALSE      66     66
#> 78        FALSE      66     66
#> 132       FALSE      66     66
#> 56        FALSE      66     66
#> 74        FALSE      66     66
#> 100       FALSE      66     66
#> 29        FALSE      66     66
#> 150       FALSE      66     66
#> 203       FALSE      66     66
#> 16        FALSE      66     66
#> 183       FALSE      66     66
#> 85        FALSE      66     66
#> 135       FALSE      66     66
#> 171       FALSE      66     66
#> 123       FALSE      66     66
#> 32        FALSE      66     66
#> 140       FALSE      66     66
#> 112       FALSE      66     66
#> 148       FALSE      66     66
#> 198       FALSE      66     66
#> 180       FALSE      66     66
#> 58        FALSE      66     66
#> 177       FALSE      66     66
Design <- oncology_scenario()[1, ]
Design
#>                                                         mu_W
#> 1 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#>                                                         mu_L
#> 1 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Sigma_W_L
#> 1 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#>                                                                        beta_prog
#> 1 0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#>                                                                   beta_switch
#> 1 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#>                                                                                     beta_death
#> 1 -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#>                                                                         beta_cens
#> 1 -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#>    k recr_interval max_duration alpha power p_trt w
#> 1 10             2            7  0.05   0.8   0.5 0
#>   aimed_for_n_per_required_event version
#> 1                              3       1
oncology_scenario_set_truevalues(oncology_scenario())
#>                                                                                                                mu_W
#> 1                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 2        0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 3  1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 4                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 5                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 6                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 7                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 8                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 9                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 10                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 11                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 12                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 13                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 14                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 15                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 16                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 17                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 18                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 19                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 20                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 21                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 22                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 23                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 24                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 25                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 26                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 27                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 28       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 29 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 30                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 31                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 32                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 33                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 34                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 35                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 36                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 37                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 38                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 39                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 40                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 41                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 42                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 43                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 44                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 45                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 46                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 47                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 48                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 49                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 50                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 51                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 52                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 53       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 54 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 55                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 56                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 57                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 58                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 59                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 60                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 61                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 62                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 63                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 64                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 65                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 66                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 67                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 68                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 69                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 70                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 71                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 72                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 73                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 74                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 75                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 76                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 77       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 78 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 79                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 80                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 81                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 82                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 83                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 84                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 85                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 86                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 87                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 88                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 89                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 90                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 91                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 92                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 93                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 94                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 95                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 96                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 97                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 98                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 99                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#>                                                                                                                mu_L
#> 1                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 2                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 3                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 4        0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 5  1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 6                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 7                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 8                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 9                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 10                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 11                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 12                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 13                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 14                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 15                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 16                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 17                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 18                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 19                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 20                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 21                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 22                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 23                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 24                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 25                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 26                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 27                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 28                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 29                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 30       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 31 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 32                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 33                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 34                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 35                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 36                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 37                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 38                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 39                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 40                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 41                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 42                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 43                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 44                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 45                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 46                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 47                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 48                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 49                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 50                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 51                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 52                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 53                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 54                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 55       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 56 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 57                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 58                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 59                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 60                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 61                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 62                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 63                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 64                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 65                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 66                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 67                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 68                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 69                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 70                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 71                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 72                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 73                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 74                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 75                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 76                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 77                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 78                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 79       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 80 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.5, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
#> 81                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 82                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 83                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 84                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 85                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 86                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 87                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 88                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 89                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 90                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 91                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 92                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 93                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 94                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 95                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 96                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 97                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 98                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#> 99                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Sigma_W_L
#> 1  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 2  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 3  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 4  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 5  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 6  1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
#> 7  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 8  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 9  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 10 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 11 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 12 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 13 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 14 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 15 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 16 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 17 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 18 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 19 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 20 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 21 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 22 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 23 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 24 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 25 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 26 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 27 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 28 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 29 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 30 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 31 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 32 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
#> 33 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 34 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 35 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 36 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 37 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 38 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 39 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 40 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 41 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 42 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 43 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 44 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 45 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 46 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 47 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 48 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 49 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 50 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 51 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 52 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 53 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 54 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 55 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 56 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 57 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
#> 58 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 59 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 60 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 61 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 62 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 63 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 64 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 65 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 66 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 67 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 68 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 69 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 70 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 71 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 72 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 73 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 74 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 75 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 76 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 77 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 78 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 79 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 80 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 81 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.6, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.7, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.9, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
#> 82 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 83 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 84 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 85 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 86 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 87 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 88 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 89 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 90 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 91 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 92 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 93 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 94 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 95 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 96 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 97 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 98 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#> 99 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0
#>                                                                          beta_prog
#> 1   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 2   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 3   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 4  0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 5  0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 6   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 7  -0.3665129, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 8    0.3266343, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 9   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, 0.0000000
#> 10   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 11 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 12  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 13  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 14  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 15  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 16  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 17  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 18 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 19  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 20  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 21  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 22  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 23  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 24  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 25  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 26  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 27  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 28  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 29  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 30 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 31 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 32  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 33 -0.3665129, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 34   0.3266343, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 35  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, 0.0000000
#> 36   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 37 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 38  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 39  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 40  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 41  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 42  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 43  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 44 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 45  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 46  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 47  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 48  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 49  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 50  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 51  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 52  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 53  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 54  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 55 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 56 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 57  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 58 -0.3665129, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 59   0.3266343, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 60  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, 0.0000000
#> 61   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 62 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 63  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 64  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 65  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 66  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 67  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 68  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 69 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 70  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 71  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 72  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 73  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 74  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 75  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 76  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 77  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 78  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 79 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 80 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 81  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 82 -0.3665129, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 83   0.3266343, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 84  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, 0.0000000
#> 85   0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 86 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 87  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 88  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 89  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 90  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 91  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 92  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 93 0.3266343, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, 0.0000000
#> 94  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 95  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 96  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 97  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 98  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#> 99  0.3266343, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000
#>                                                                    beta_switch
#> 1  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 2  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 3  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 4  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 5  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 6  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 7  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 8  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 9  0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 10 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 11 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 12 2.1972246, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 13 0.0000000, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 14 0.0000000, 0.4054651, 0.4054651, 2.3483475, 0.0000000, 0.0000000, 0.0000000
#> 15 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 16 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 17 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 18 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 19 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 20 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 21 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 22 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 23 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 24 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 25 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 26 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 27 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 28 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 29 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 30 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 31 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 32 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 33 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 34 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 35 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 36 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 37 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 38 2.1972246, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 39 0.0000000, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 40 0.0000000, 0.4054651, 0.4054651, 2.3483475, 0.0000000, 0.0000000, 0.0000000
#> 41 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 42 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 43 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 44 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 45 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 46 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 47 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 48 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 49 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 50 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 51 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 52 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 53 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 54 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 55 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 56 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 57 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 58 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 59 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 60 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 61 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 62 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 63 2.1972246, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 64 0.0000000, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 65 0.0000000, 0.4054651, 0.4054651, 2.3483475, 0.0000000, 0.0000000, 0.0000000
#> 66 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 67 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 68 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 69 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 70 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 71 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 72 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 73 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 74 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 75 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 76 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 77 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 78 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 79 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 80 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 81 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 82 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 83 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 84 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 85 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 86 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 87 2.1972246, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 88 0.0000000, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 89 0.0000000, 0.4054651, 0.4054651, 2.3483475, 0.0000000, 0.0000000, 0.0000000
#> 90 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 91 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 92 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 93 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.4054651, 0.0000000, 0.0000000
#> 94 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 95 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 96 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 97 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 98 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 99 0.0000000, 0.4054651, 0.4054651, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#>                                                                                       beta_death
#> 1   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 2   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 3   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 4  -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, -0.6931472, -0.6931472
#> 5  -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, -0.6931472, -0.6931472
#> 6   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 7   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 8   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 9   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 10  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 11 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, -0.6931472, -0.6931472
#> 12  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 13  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 14  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 15 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, -0.6931472, -0.6931472
#> 16  -1.0596601, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 17   -1.7528073, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 18 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.6931472, -0.6931472, -0.6931472
#> 19  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.2876821, -0.6931472
#> 20   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, 0.0000000, -0.6931472
#> 21  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 22  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 23  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 24  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 25  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 26  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.6931472, -0.6931472, -0.6931472
#> 27  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 28  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 29  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 30 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.2876821, -0.2876821, -0.2876821
#> 31 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.2876821, -0.2876821, -0.2876821
#> 32  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 33  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 34  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 35  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 36  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 37 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.2876821, -0.2876821, -0.2876821
#> 38  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 39  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 40  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 41 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.2876821, -0.2876821, -0.2876821
#> 42  -1.0596601, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 43   -1.7528073, -0.6931472, 0.0000000, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 44 -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, -0.2876821, -0.2876821, -0.2876821
#> 45   -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, 0.0000000, -0.2876821
#> 46  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 47  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 48  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 49  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 50  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 51  -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, -0.2876821, -0.2876821, -0.2876821
#> 52    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 53    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 54    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 55   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.6931472
#> 56   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.6931472
#> 57    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 58    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 59    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 60    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 61    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 62   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.6931472
#> 63    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 64    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 65    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 66   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.6931472
#> 67    -1.0596601, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 68     -1.7528073, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 69   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.6931472
#> 70    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 71    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 72    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 73    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 74    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 75    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.6931472
#> 76    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 77    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 78    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 79   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.2876821
#> 80   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.2876821
#> 81    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 82    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 83    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 84    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 85    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 86   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.2876821
#> 87    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 88    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 89    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 90   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.2876821
#> 91    -1.0596601, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 92     -1.7528073, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 93   -1.7528073, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, -0.2876821
#> 94    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 95    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 96    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 97    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 98    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#> 99    -1.7528073, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -0.2876821
#>                                                                                     beta_cens
#> 1             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 2             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 3             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 4             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 5             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 6             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 7             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 8             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 9             -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 10            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 11            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 12            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 13            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 14            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 15            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 16            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 17            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 18            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 19            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 20            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 21                                                                  -Inf, 0, 0, 0, 0, 0, 0, 0
#> 22   -3.6762473, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 23  -3.6762473, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 24 -3.6762473, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, 0.0000000
#> 25  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 26  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000
#> 27            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 28            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 29            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 30            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 31            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 32            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 33            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 34            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 35            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 36            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 37            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 38            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 39            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 40            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 41            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 42            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 43            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 44            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 45            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 46                                                                  -Inf, 0, 0, 0, 0, 0, 0, 0
#> 47   -3.6762473, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 48  -3.6762473, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 49 -3.6762473, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, 0.0000000
#> 50  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 51  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000
#> 52            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 53            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 54            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 55            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 56            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 57            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 58            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 59            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 60            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 61            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 62            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 63            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 64            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 65            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 66            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 67            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 68            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 69            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 70                                                                  -Inf, 0, 0, 0, 0, 0, 0, 0
#> 71   -3.6762473, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 72  -3.6762473, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 73 -3.6762473, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, 0.0000000
#> 74  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 75  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000
#> 76            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 77            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 78            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 79            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 80            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 81            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 82            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 83            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 84            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 85            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 86            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 87            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 88            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 89            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 90            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 91            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 92            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 93            -3.676247, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
#> 94                                                                  -Inf, 0, 0, 0, 0, 0, 0, 0
#> 95   -3.6762473, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 96  -3.6762473, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 97 -3.6762473, -0.6931472, -0.6931472, 0.0000000, -0.6931472, 0.0000000, 0.0000000, 0.0000000
#> 98  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000
#> 99  -2.2503673, -0.6931472, -0.6931472, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 1.0000000
#>     k recr_interval max_duration alpha power p_trt w
#> 1  10             2            7  0.05   0.8   0.5 0
#> 2  10             2            7  0.05   0.8   0.5 0
#> 3  10             2            7  0.05   0.8   0.5 0
#> 4  10             2            7  0.05   0.8   0.5 0
#> 5  10             2            7  0.05   0.8   0.5 0
#> 6  10             2            7  0.05   0.8   0.5 0
#> 7  10             2            7  0.05   0.8   0.5 0
#> 8  10             2            7  0.05   0.8   0.5 0
#> 9  10             2            7  0.05   0.8   0.5 0
#> 10 10             2            7  0.05   0.8   0.5 0
#> 11 10             2            7  0.05   0.8   0.5 0
#> 12 10             2            7  0.05   0.8   0.5 0
#> 13 10             2            7  0.05   0.8   0.5 0
#> 14 10             2            7  0.05   0.8   0.5 0
#> 15 10             2            7  0.05   0.8   0.5 0
#> 16 10             2            7  0.05   0.8   0.5 0
#> 17 10             2            7  0.05   0.8   0.5 0
#> 18 10             2            7  0.05   0.8   0.5 0
#> 19 10             2            7  0.05   0.8   0.5 0
#> 20 10             2            7  0.05   0.8   0.5 0
#> 21 10             2            7  0.05   0.8   0.5 0
#> 22 10             2            7  0.05   0.8   0.5 0
#> 23 10             2            7  0.05   0.8   0.5 0
#> 24 10             2            7  0.05   0.8   0.5 0
#> 25 10             2            7  0.05   0.8   0.5 0
#> 26 10             2            7  0.05   0.8   0.5 0
#> 27 10             2            7  0.05   0.8   0.5 0
#> 28 10             2            7  0.05   0.8   0.5 0
#> 29 10             2            7  0.05   0.8   0.5 0
#> 30 10             2            7  0.05   0.8   0.5 0
#> 31 10             2            7  0.05   0.8   0.5 0
#> 32 10             2            7  0.05   0.8   0.5 0
#> 33 10             2            7  0.05   0.8   0.5 0
#> 34 10             2            7  0.05   0.8   0.5 0
#> 35 10             2            7  0.05   0.8   0.5 0
#> 36 10             2            7  0.05   0.8   0.5 0
#> 37 10             2            7  0.05   0.8   0.5 0
#> 38 10             2            7  0.05   0.8   0.5 0
#> 39 10             2            7  0.05   0.8   0.5 0
#> 40 10             2            7  0.05   0.8   0.5 0
#> 41 10             2            7  0.05   0.8   0.5 0
#> 42 10             2            7  0.05   0.8   0.5 0
#> 43 10             2            7  0.05   0.8   0.5 0
#> 44 10             2            7  0.05   0.8   0.5 0
#> 45 10             2            7  0.05   0.8   0.5 0
#> 46 10             2            7  0.05   0.8   0.5 0
#> 47 10             2            7  0.05   0.8   0.5 0
#> 48 10             2            7  0.05   0.8   0.5 0
#> 49 10             2            7  0.05   0.8   0.5 0
#> 50 10             2            7  0.05   0.8   0.5 0
#> 51 10             2            7  0.05   0.8   0.5 0
#> 52 10             2            7  0.05   0.8   0.5 0
#> 53 10             2            7  0.05   0.8   0.5 0
#> 54 10             2            7  0.05   0.8   0.5 0
#> 55 10             2            7  0.05   0.8   0.5 0
#> 56 10             2            7  0.05   0.8   0.5 0
#> 57 10             2            7  0.05   0.8   0.5 0
#> 58 10             2            7  0.05   0.8   0.5 0
#> 59 10             2            7  0.05   0.8   0.5 0
#> 60 10             2            7  0.05   0.8   0.5 0
#> 61 10             2            7  0.05   0.8   0.5 0
#> 62 10             2            7  0.05   0.8   0.5 0
#> 63 10             2            7  0.05   0.8   0.5 0
#> 64 10             2            7  0.05   0.8   0.5 0
#> 65 10             2            7  0.05   0.8   0.5 0
#> 66 10             2            7  0.05   0.8   0.5 0
#> 67 10             2            7  0.05   0.8   0.5 0
#> 68 10             2            7  0.05   0.8   0.5 0
#> 69 10             2            7  0.05   0.8   0.5 0
#> 70 10             2            7  0.05   0.8   0.5 0
#> 71 10             2            7  0.05   0.8   0.5 0
#> 72 10             2            7  0.05   0.8   0.5 0
#> 73 10             2            7  0.05   0.8   0.5 0
#> 74 10             2            7  0.05   0.8   0.5 0
#> 75 10             2            7  0.05   0.8   0.5 0
#> 76 10             2            7  0.05   0.8   0.5 0
#> 77 10             2            7  0.05   0.8   0.5 0
#> 78 10             2            7  0.05   0.8   0.5 0
#> 79 10             2            7  0.05   0.8   0.5 0
#> 80 10             2            7  0.05   0.8   0.5 0
#> 81 10             2            7  0.05   0.8   0.5 0
#> 82 10             2            7  0.05   0.8   0.5 0
#> 83 10             2            7  0.05   0.8   0.5 0
#> 84 10             2            7  0.05   0.8   0.5 0
#> 85 10             2            7  0.05   0.8   0.5 0
#> 86 10             2            7  0.05   0.8   0.5 0
#> 87 10             2            7  0.05   0.8   0.5 0
#> 88 10             2            7  0.05   0.8   0.5 0
#> 89 10             2            7  0.05   0.8   0.5 0
#> 90 10             2            7  0.05   0.8   0.5 0
#> 91 10             2            7  0.05   0.8   0.5 0
#> 92 10             2            7  0.05   0.8   0.5 0
#> 93 10             2            7  0.05   0.8   0.5 0
#> 94 10             2            7  0.05   0.8   0.5 0
#> 95 10             2            7  0.05   0.8   0.5 0
#> 96 10             2            7  0.05   0.8   0.5 0
#> 97 10             2            7  0.05   0.8   0.5 0
#> 98 10             2            7  0.05   0.8   0.5 0
#> 99 10             2            7  0.05   0.8   0.5 0
#>    aimed_for_n_per_required_event version ev_soll
#> 1                               3       1      66
#> 2                               3       1      66
#> 3                               3       1      66
#> 4                               3       1      66
#> 5                               3       1      66
#> 6                               3       1      66
#> 7                               3       1      66
#> 8                               3       1      66
#> 9                               3       1      66
#> 10                              3       1      66
#> 11                              3       1      66
#> 12                              3       1      66
#> 13                              3       1      66
#> 14                              3       1      66
#> 15                              3       1      66
#> 16                              3       1      66
#> 17                              3       1      66
#> 18                              3       1      66
#> 19                              3       1      66
#> 20                              3       1      66
#> 21                              3       1      66
#> 22                              3       1      66
#> 23                              3       1      66
#> 24                              3       1      66
#> 25                              3       1      66
#> 26                              3       1      66
#> 27                              3       1     380
#> 28                              3       1     380
#> 29                              3       1     380
#> 30                              3       1     380
#> 31                              3       1     380
#> 32                              3       1     380
#> 33                              3       1     380
#> 34                              3       1     380
#> 35                              3       1     380
#> 36                              3       1     380
#> 37                              3       1     380
#> 38                              3       1     380
#> 39                              3       1     380
#> 40                              3       1     380
#> 41                              3       1     380
#> 42                              3       1     380
#> 43                              3       1     380
#> 44                              3       1     380
#> 45                              3       1     380
#> 46                              3       1     380
#> 47                              3       1     380
#> 48                              3       1     380
#> 49                              3       1     380
#> 50                              3       1     380
#> 51                              3       1     380
#> 52                              3       1      66
#> 53                              3       1      66
#> 54                              3       1      66
#> 55                              3       1      66
#> 56                              3       1      66
#> 57                              3       1      66
#> 58                              3       1      66
#> 59                              3       1      66
#> 60                              3       1      66
#> 61                              3       1      66
#> 62                              3       1      66
#> 63                              3       1      66
#> 64                              3       1      66
#> 65                              3       1      66
#> 66                              3       1      66
#> 67                              3       1      66
#> 68                              3       1      66
#> 69                              3       1      66
#> 70                              3       1      66
#> 71                              3       1      66
#> 72                              3       1      66
#> 73                              3       1      66
#> 74                              3       1      66
#> 75                              3       1      66
#> 76                              3       1     380
#> 77                              3       1     380
#> 78                              3       1     380
#> 79                              3       1     380
#> 80                              3       1     380
#> 81                              3       1     380
#> 82                              3       1     380
#> 83                              3       1     380
#> 84                              3       1     380
#> 85                              3       1     380
#> 86                              3       1     380
#> 87                              3       1     380
#> 88                              3       1     380
#> 89                              3       1     380
#> 90                              3       1     380
#> 91                              3       1     380
#> 92                              3       1     380
#> 93                              3       1     380
#> 94                              3       1     380
#> 95                              3       1     380
#> 96                              3       1     380
#> 97                              3       1     380
#> 98                              3       1     380
#> 99                              3       1     380
```
