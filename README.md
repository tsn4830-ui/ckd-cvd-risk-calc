# 腎臟與心血管風險計算器

繁體中文的臨床風險計算工具，全下拉選單操作，預設台灣常用單位（mg/dL、mg/g、HbA1c %）。

| 模型 | 用途 | 適用對象 |
|---|---|---|
| **KFRE** | 2 年／5 年腎衰竭（需腎臟替代治療）風險 | CKD G3–G5 |
| **SCORE2 + CKD Add-on** | 10 年心血管風險 | 40–69 歲，無心血管疾病 |
| **SCORE2-OP + CKD Add-on** | 10 年心血管風險 | 70–89 歲，無心血管疾病 |
| **SCORE2-Diabetes** | 10 年心血管風險 | 第 2 型糖尿病 40–69 歲，無心血管疾病 |

純靜態網站（HTML/CSS/JS），所有計算在瀏覽器本機完成，不傳送任何資料。

## 演算法驗證

- **KFRE**：係數與基準存活率與官方計算器 kidneyfailurerisk.com 的程式碼逐項比對，4 變數與 8 變數、北美／非北美校正、2 年／5 年全部吻合。
- **SCORE2 / SCORE2-OP + CKD Add-on**：與 CKD-PC 官方計算器（ckdpcrisk.org）針對 40–85 歲、雙性別、四風險地區、eGFR 15–120、uACR 5–1500 共 **280 組數值逐一比對，完全吻合**。
- **SCORE2-Diabetes**：以原論文內文的四組範例值驗證吻合。

## 文獻

- Tangri N, et al. *JAMA* 2011;305:1553-9；*JAMA* 2016;315:164-74
- SCORE2 working group. *Eur Heart J* 2021;42:2439-54
- SCORE2-OP working group. *Eur Heart J* 2021;42:2455-67
- Matsushita K, et al. *Eur J Prev Cardiol* 2023;30:8-16
- SCORE2-Diabetes working group. *Eur Heart J* 2023;44:2544-56
- KDIGO 2024 CKD Guideline；ESC 2021 CVD Prevention；ESC 2023 Diabetes

## 免責

僅供醫療專業人員參考，不能取代臨床判斷。SCORE2 系列以歐洲族群校正，台灣不屬於其四個風險地區，解讀時請留意。
