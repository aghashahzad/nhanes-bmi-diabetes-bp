# BMI, diabetes and high blood pressure in US adults (NHANES 2021-2023)

A cross-sectional analysis of public survey data using logistic regression, written in Python.

## Question
Are overweight and obesity associated with diabetes and with high measured blood pressure in adults, after adjusting for age and sex?

## Data
- Source: CDC National Center for Health Statistics, National Health and Nutrition Examination Survey (NHANES), August 2021-August 2023.
- Download page: https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx
- Files used (all from the 2021-2023 cycle): `DEMO_L` (demographics), `BMX_L` (body measures), `GHB_L` (HbA1c), `DIQ_L` (diabetes questionnaire), `BPXO_L` (blood pressure).
- To rerun: download the five .xpt files from the page above into the same folder as `analysis.py`.

## Definitions
- **Diabetes:** told by a doctor that they have diabetes, or HbA1c of 6.5% or higher.
- **High measured blood pressure:** mean systolic 140 mmHg or higher, or mean diastolic 90 mmHg or higher, using the available readings from one examination. This is not a clinical diagnosis of hypertension, because medication use was not included.
- **BMI categories:** normal 18.5-24.9 (reference), overweight 25-29.9, obesity 30 or higher. Underweight participants (BMI under 18.5) were excluded.

## Sample
| Step | n |
|---|---|
| All participants | 11,933 |
| Age 20 or older | 7,809 |
| Not pregnant | 7,768 |
| BMI measured | 5,929 |
| Not underweight | 5,847 |

Participants with a missing outcome were left out of that model: 5,624 for diabetes and 5,689 for blood pressure.

## Method
1. Merged the five files by participant ID.
2. Calculated the share of participants with each outcome in each BMI category.
3. Fitted logistic regression models: crude, adjusted for age and sex, and BMI per 5 units. The regression is written out step by step in `analysis.py` and checked against scikit-learn (coefficients agree to within 0.0001).

## Results

**Prevalence by BMI category (unweighted sample percentages)**

| Outcome | Normal | Overweight | Obesity |
|---|---|---|---|
| Diabetes | 8.0% | 14.4% | 26.0% |
| High measured blood pressure | 16.7% | 18.5% | 20.4% |

**Odds ratios (95% confidence interval), reference: normal weight**

| Outcome | Category | Crude OR | Adjusted OR (age, sex) |
|---|---|---|---|
| Diabetes | Overweight | 1.92 (1.52-2.42) | 1.81 (1.43-2.30) |
| Diabetes | Obesity | 4.03 (3.26-4.98) | 4.30 (3.46-5.35) |
| High blood pressure | Overweight | 1.13 (0.94-1.35) | 1.05 (0.87-1.26) |
| High blood pressure | Obesity | 1.28 (1.08-1.51) | 1.25 (1.05-1.49) |

Per 5 units of BMI, the adjusted OR was 1.49 (1.42-1.57) for diabetes and 1.11 (1.05-1.16) for high blood pressure.

![Prevalence](chart1_prevalence_by_bmi.png)
![Odds ratios](chart2_odds_ratios.png)

## Findings
- Diabetes was much more common with higher BMI: 8.0% in normal weight versus 26.0% with obesity. Adults with obesity had about 4.3 times the odds of diabetes compared with normal weight, after adjusting for age and sex.
- The association with high measured blood pressure was weaker. The adjusted odds ratio was 1.25 for obesity, and the interval for overweight included 1.
- Older age was also strongly associated with both outcomes.

## Limitations
- **Cross-sectional design:** it cannot show that higher BMI causes either outcome, or the order in which they developed.
- **Unweighted analysis:** NHANES uses a complex sample design. I did not apply survey weights, so the percentages describe this sample and are not national estimates.
- **Limited adjustment:** only age and sex were included. Race and ethnicity, income, smoking, diet and physical activity could change the results.
- **Outcome definitions:** blood pressure came from one examination without medication information, and diabetes does not distinguish type 1 from type 2.
- **Complete cases only,** and BMI is an imperfect measure of body fat.
- US data, so the results may not apply to other populations.

## Files
- `analysis.py`: full script
- `sample_flow.csv`, `prevalence_by_bmi.csv`, `odds_ratios.csv`: result tables
- `chart1_prevalence_by_bmi.png`, `chart2_odds_ratios.png`: charts

## Author
Agha Shahzad Hussain, BS Radiologic Technology (LUMHS, Pakistan, 2025). Self-directed project in epidemiology and data analysis, in preparation for postgraduate study in public health.
