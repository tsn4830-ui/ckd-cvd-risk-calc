# 腎臟與心血管風險計算器

繁體中文的臨床風險計算工具，全下拉選單操作，預設台灣常用單位（mg/dL、mg/g、HbA1c %）。

| 模型 | 用途 | 適用對象 |
|---|---|---|
| **KFRE** | 2 年／5 年腎衰竭（需腎臟替代治療）風險 | CKD G3–G5 |
| **SCORE2 + CKD Add-on** | 10 年心血管風險 | 40–69 歲，無心血管疾病 |
| **SCORE2-OP + CKD Add-on** | 10 年心血管風險 | 70–89 歲，無心血管疾病 |
| **SCORE2-Diabetes** | 10 年心血管風險 | 第 2 型糖尿病 40–69 歲，無心血管疾病 |
| **SCORE2-Asia-Pacific** | 10 年心血管風險（亞太校正） | 40–69 歲，無心血管疾病且無糖尿病 |

純靜態網站（HTML/CSS/JS），所有計算在瀏覽器本機完成，不傳送任何資料。

其他功能：
- **KDIGO CKD 分級熱區圖**：KFRE 分頁會依 eGFR／uACR 把病人標在 KDIGO 2024 Figure 13 的 6×3 格子上，格內數字為該格建議的每年監測次數。
- **列印評估摘要**：每個分頁都有「列印此次評估」，輸出 A4 單頁（輸入條件、結果、判讀建議、熱區圖、姓名病歷號空白欄），可夾病歷或給病人帶走。
- **單位切換**：膽固醇 mg/dL ↔ mmol/L、HbA1c % ↔ mmol/mol、uACR mg/g ↔ mg/mmol、血清白蛋白 g/dL ↔ g/L。
- **CKD-EPI 2021 eGFR 換算小工具**。

## 演算法驗證

- **KFRE**：係數與基準存活率與官方計算器 kidneyfailurerisk.com 的程式碼逐項比對，4 變數與 8 變數、北美／非北美校正、2 年／5 年全部吻合。
- **SCORE2 / SCORE2-OP + CKD Add-on**：與 CKD-PC 官方計算器（ckdpcrisk.org）針對 40–85 歲、雙性別、四風險地區、eGFR 15–120、uACR 5–1500 共 **280 組數值逐一比對，完全吻合**。
- **SCORE2-Diabetes**：以原論文內文的四組範例值驗證吻合。
- **SCORE2-Asia-Pacific**：校正因子取自原論文 Supplementary Table 5（Cambridge repository 的 accepted 版），並以內文範例（50 歲、SBP 140、TC 5.5、HDL 1.3 之男女，四區）驗證吻合。

### 台灣屬於哪個亞太風險區？

原論文的國家清單（Supplementary Table 4）**沒有台灣**——分區依 WHO Global Health Estimates，而 WHO 統計不含台灣。分區門檻為年齡性別標準化心血管死亡率每十萬人年：低 <100、中 100–<150、高 150–<300、極高 ≥300。以同一資料源（IHME GBD）比較：日本 77、新加坡 85、**台灣 90**、泰國 103（前兩者論文列為低風險區，泰國為中風險區）。因此本站將台灣對應到**低風險區**，但同時顯示低與中風險區作為判讀區間。此為本站依公開資料所做的對應，非原論文的官方分類。

## 文獻

- Tangri N, et al. *JAMA* 2011;305:1553-9；*JAMA* 2016;315:164-74
- SCORE2 working group. *Eur Heart J* 2021;42:2439-54
- SCORE2-OP working group. *Eur Heart J* 2021;42:2455-67
- Matsushita K, et al. *Eur J Prev Cardiol* 2023;30:8-16
- SCORE2-Diabetes working group. *Eur Heart J* 2023;44:2544-56
- SCORE2 Asia-Pacific collaborators. *Eur Heart J* 2025;46:702-15
- KDIGO 2024 CKD Guideline；ESC 2021 CVD Prevention；ESC 2023 Diabetes

## 免責

僅供醫療專業人員參考，不能取代臨床判斷。SCORE2／SCORE2-OP／SCORE2-Diabetes 以歐洲族群校正；亞太族群請優先參考 SCORE2-Asia-Pacific 分頁，並留意台灣的區域歸屬是推估而非官方分類。
