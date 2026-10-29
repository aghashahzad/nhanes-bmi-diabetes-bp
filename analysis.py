"""
BMI, diabetes and high blood pressure in US adults (NHANES August 2021-August 2023)

Data: CDC National Center for Health Statistics, NHANES 2021-2023 public files
      https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx
Files (place the .xpt files in the same folder as this script):
      DEMO_L (demographics), BMX_L (body measures), GHB_L (HbA1c),
      DIQ_L (diabetes questionnaire), BPXO_L (blood pressure)

Method: logistic regression, written out with numpy (iteratively reweighted
least squares) so each step is visible, and checked against scikit-learn.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# 1. Load and merge the five files on SEQN (participant ID)
# ---------------------------------------------------------------------
files = ["DEMO_L", "BMX_L", "GHB_L", "DIQ_L", "BPXO_L"]
data = {f: pd.read_sas(f + ".xpt", format="xport") for f in files}
for f, d in data.items():
    print(f, d.shape)

df = data["DEMO_L"][["SEQN", "RIDAGEYR", "RIAGENDR", "RIDEXPRG"]]
for f, cols in [("BMX_L", ["SEQN", "BMXBMI"]),
                ("GHB_L", ["SEQN", "LBXGH"]),
                ("DIQ_L", ["SEQN", "DIQ010"]),
                ("BPXO_L", ["SEQN", "BPXOSY1", "BPXOSY2", "BPXOSY3",
                            "BPXODI1", "BPXODI2", "BPXODI3"])]:
    df = df.merge(data[f][cols], on="SEQN", how="left")
print("Merged participants:", len(df))

# ---------------------------------------------------------------------
# 2. Define variables
# ---------------------------------------------------------------------
df = df.rename(columns={"RIDAGEYR": "age", "BMXBMI": "bmi", "LBXGH": "hba1c"})
df["female"] = (df["RIAGENDR"] == 2).astype(int)

# Diabetes: doctor-told diabetes OR HbA1c >= 6.5%
told = df["DIQ010"] == 1
no_told = df["DIQ010"].isin([2, 3])
df["diabetes"] = np.nan
df.loc[no_told & (df["hba1c"] < 6.5), "diabetes"] = 0
df.loc[told | (df["hba1c"] >= 6.5), "diabetes"] = 1

# High measured blood pressure: mean of available readings, SBP >= 140 or DBP >= 90
sbp = df[["BPXOSY1", "BPXOSY2", "BPXOSY3"]].mean(axis=1)
dbp = df[["BPXODI1", "BPXODI2", "BPXODI3"]].mean(axis=1)
df["high_bp"] = np.where(sbp.isna() | dbp.isna(), np.nan, ((sbp >= 140) | (dbp >= 90)).astype(float))

# BMI category (WHO cut-offs)
df["bmi_cat"] = pd.cut(df["bmi"], [0, 18.5, 25, 30, 100], right=False,
                       labels=["Underweight", "Normal", "Overweight", "Obesity"])

# ---------------------------------------------------------------------
# 3. Select the analysis sample and keep a record of every exclusion
# ---------------------------------------------------------------------
flow = [("All participants", len(df))]
a = df[df["age"] >= 20];                      flow.append(("Age 20 years or older", len(a)))
a = a[a["RIDEXPRG"] != 1];                    flow.append(("Not pregnant", len(a)))
a = a[a["bmi"].notna()];                      flow.append(("BMI measured", len(a)))
a = a[a["bmi_cat"] != "Underweight"];         flow.append(("Not underweight (BMI >= 18.5)", len(a)))
a = a.copy()
a["bmi_cat"] = a["bmi_cat"].cat.remove_unused_categories()
pd.DataFrame(flow, columns=["step", "n"]).to_csv("sample_flow.csv", index=False)
print(pd.DataFrame(flow, columns=["step", "n"]))

# ---------------------------------------------------------------------
# 4. Prevalence by BMI category
# ---------------------------------------------------------------------
prev = []
for outcome, label in [("diabetes", "Diabetes"), ("high_bp", "High measured blood pressure")]:
    sub = a[a[outcome].notna()]
    t = sub.groupby("bmi_cat", observed=True)[outcome].agg(["size", "mean"])
    t["mean"] = (t["mean"] * 100).round(1)
    t = t.reset_index()
    t.insert(0, "outcome", label)
    t.columns = ["outcome", "bmi_category", "n", "prevalence_percent"]
    prev.append(t)
prev = pd.concat(prev)
print(prev)
prev.to_csv("prevalence_by_bmi.csv", index=False)

# ---------------------------------------------------------------------
# 5. Logistic regression (IRLS, written out so each step is visible)
# ---------------------------------------------------------------------
def logit(X, y, n_iter=50):
    """Return coefficients and standard errors for a logistic regression."""
    X = np.asarray(X, dtype=float); y = np.asarray(y, dtype=float)
    beta = np.zeros(X.shape[1])
    for _ in range(n_iter):
        p = 1 / (1 + np.exp(-X @ beta))
        W = p * (1 - p)
        grad = X.T @ (y - p)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, grad)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            break
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    return beta, np.sqrt(np.diag(cov))

def odds_ratios(sub, outcome, predictors, label):
    X = sub[predictors].astype(float).copy()
    X.insert(0, "const", 1.0)
    beta, se = logit(X, sub[outcome])
    out = pd.DataFrame({
        "outcome": label, "term": X.columns, "OR": np.exp(beta),
        "CI_low": np.exp(beta - 1.96 * se), "CI_high": np.exp(beta + 1.96 * se),
    })
    return out[out["term"] != "const"], X, sub[outcome].values, beta

a["overweight"] = (a["bmi_cat"] == "Overweight").astype(int)
a["obesity"] = (a["bmi_cat"] == "Obesity").astype(int)
a["bmi_per5"] = a["bmi"] / 5
a["age_per10"] = a["age"] / 10

results = []
for outcome, label in [("diabetes", "Diabetes"), ("high_bp", "High measured blood pressure")]:
    sub = a[a[outcome].notna()]
    n = len(sub)
    for model, preds in [("Crude", ["overweight", "obesity"]),
                         ("Adjusted for age and sex", ["overweight", "obesity", "age_per10", "female"])]:
        r, X, y, beta = odds_ratios(sub, outcome, preds, label)
        r["model"] = model; r["n"] = n
        results.append(r)
    r, X, y, beta = odds_ratios(sub, outcome, ["bmi_per5", "age_per10", "female"], label)
    r["model"] = "Adjusted, BMI per 5 units"; r["n"] = n
    results.append(r)

res = pd.concat(results)[["outcome", "model", "term", "n", "OR", "CI_low", "CI_high"]]
res[["OR", "CI_low", "CI_high"]] = res[["OR", "CI_low", "CI_high"]].round(2)
print(res.to_string(index=False))
res.to_csv("odds_ratios.csv", index=False)

# Check against scikit-learn (no penalty) on one model
from sklearn.linear_model import LogisticRegression
sub = a[a["diabetes"].notna()]
Xs = sub[["overweight", "obesity", "age_per10", "female"]].astype(float)
sk = LogisticRegression(C=np.inf, max_iter=1000).fit(Xs, sub["diabetes"])
mine = odds_ratios(sub, "diabetes", ["overweight", "obesity", "age_per10", "female"], "check")[0]
print("Check vs scikit-learn (max difference in log-OR):",
      float(np.max(np.abs(np.log(mine["OR"].values) - sk.coef_[0]))))

# ---------------------------------------------------------------------
# 6. Charts
# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)
for ax, (label, color) in zip(axes, [("Diabetes", "#1F3864"), ("High measured blood pressure", "#C55A11")]):
    t = prev[prev["outcome"] == label]
    bars = ax.bar(t["bmi_category"], t["prevalence_percent"], color=color)
    for b, v in zip(bars, t["prevalence_percent"]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v}%", ha="center")
    ax.set_ylim(0, t["prevalence_percent"].max() * 1.18)
    ax.set_title(label)
    ax.set_ylabel("Prevalence (%)")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Prevalence by BMI category, US adults 20+ (NHANES 2021-2023)")
fig.text(0.01, 0.01, "Source: CDC NHANES August 2021-August 2023. Unweighted sample percentages.",
         fontsize=7, color="gray")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig("chart1_prevalence_by_bmi.png", dpi=150)

adj = res[(res["model"] == "Adjusted for age and sex") & (res["term"].isin(["overweight", "obesity"]))].copy()
adj["label"] = adj["outcome"] + " - " + adj["term"].map({"overweight": "Overweight", "obesity": "Obesity"})
adj = adj.iloc[::-1]
fig, ax = plt.subplots(figsize=(8, 4))
ypos = np.arange(len(adj))
ax.errorbar(adj["OR"], ypos, xerr=[adj["OR"] - adj["CI_low"], adj["CI_high"] - adj["OR"]],
            fmt="o", color="#1F3864", capsize=4)
ax.axvline(1, color="gray", linestyle="--")
ax.set_yticks(ypos); ax.set_yticklabels(adj["label"])
ax.set_xlabel("Adjusted odds ratio (95% CI), reference: normal weight")
ax.set_title("Adjusted odds ratios by BMI category")
for y_, o_, hi_ in zip(ypos, adj["OR"], adj["CI_high"]):
    ax.text(hi_ + 0.08, y_, f"OR {o_:.2f}", ha="left", va="center", fontsize=9)
ax.set_xlim(0.5, 6.2)
ax.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.01, "Source: CDC NHANES August 2021-August 2023. Adjusted for age and sex.",
         fontsize=7, color="gray")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig("chart2_odds_ratios.png", dpi=150)
print("Done.")
