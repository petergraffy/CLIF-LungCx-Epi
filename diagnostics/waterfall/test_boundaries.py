import pandas as pd
from run_waterfall import bounded_waterfall

t0 = pd.Timestamp("2020-01-01",tz="UTC")
raw = pd.DataFrame([
    {"hospitalization_id":"a","recorded_dttm":t0,"device_category":"IMV","device_name":"Vent","fio2_set":.4},
    {"hospitalization_id":"a","recorded_dttm":t0+pd.Timedelta(hours=4),"device_category":"Room Air","device_name":"Room Air","fio2_set":.21},
    {"hospitalization_id":"b","recorded_dttm":t0,"lpm_set":2.0},
    {"hospitalization_id":"c","recorded_dttm":t0,"lpm_set":2.0},
    {"hospitalization_id":"c","recorded_dttm":t0+pd.Timedelta(hours=2),"device_category":"NIPPV","device_name":"BiPAP","fio2_set":.6},
    {"hospitalization_id":"d","recorded_dttm":t0,"device_category":"IMV","device_name":"Vent","fio2_set":.5},
    {"hospitalization_id":"d","recorded_dttm":t0+pd.Timedelta(hours=2),"device_category":"Room Air","device_name":"Room Air","fio2_set":.21},
])
bounds = pd.DataFrame({"hospitalization_id":list("abcde"),"start":[t0]*5,
                       "end":[t0+pd.Timedelta(hours=x) for x in (3,3,3,1.5,3)]})
out, _, _ = bounded_waterfall(raw,bounds)
assert set(out.hospitalization_id)==set("abcde")
assert out.loc[out.hospitalization_id=="a","device_category"].eq("imv").all()
assert out.loc[out.hospitalization_id=="b","device_category"].isna().all()
assert out.loc[(out.hospitalization_id=="c") & (out.recorded_dttm<t0+pd.Timedelta(hours=2)),"device_category"].isna().all()
assert out.loc[(out.hospitalization_id=="c") & (out.recorded_dttm<t0+pd.Timedelta(hours=2)),"fio2_set"].isna().all()
assert out.loc[out.hospitalization_id=="d","recorded_dttm"].lt(t0+pd.Timedelta(hours=1.5)).all()
assert out.loc[out.hospitalization_id=="d","device_category"].eq("imv").all()
assert out.loc[out.hospitalization_id=="e","device_category"].isna().all()
assert out.loc[(out.hospitalization_id=="b") & (out.recorded_dttm==t0),"lpm_set"].eq(2).all()
assert out.loc[out.is_scaffold.fillna(False),"raw_device_observed_at_timestamp"].eq(False).all()
print("PASS: exact landmark and terminal boundaries, no numeric backfill or unsupported device inference, retention of unobserved encounters and real boundary rows")
