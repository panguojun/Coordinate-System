using UnityEngine;
using System;

public class Life : MonoBehaviour
{
    private int lifespan; // 寿命
    private string status; // 状态（活着/死亡）
    private long boundaryStart; // 有效时间段的起始时间（毫秒）
    private long boundaryEnd; // 有效时间段的结束时间（毫秒）

    public Life(int lifespan, long boundaryStart, long boundaryEnd)
    {
        this.lifespan = lifespan;
        this.boundaryStart = boundaryStart;
        this.boundaryEnd = boundaryEnd;
        this.status = "活着";
    }

    public void UpdateStatus(long currentDate)
    {
         if (currentTime < boundaryStart || currentTime > boundaryEnd)
        {
            this.status = "死亡";
        }
        else if ((currentTime - boundaryStart) >= lifespan)
        {
            this.status = "死亡";
        }
        else
        {
            this.status = "活着";
        }
    }

    public string GetStatus()
    {
        return this.status;
    }

    // 重载加法运算符
    public static Life operator +(Life life, bool value)
    {
        if (value)
        {
            life.boundaryStart = life.boundaryStart.AddDays(1);
            life.boundaryEnd = life.boundaryEnd.AddDays(1);
        }
        else
        {
            life.boundaryStart = life.boundaryStart.AddDays(-1);
            life.boundaryEnd = life.boundaryEnd.AddDays(-1);
        }

        return life;
    }

    // 重载减法运算符
    public static Life operator -(Life life, bool value)
    {
        if (value)
        {
            life.boundaryStart = life.boundaryStart.AddDays(-1);
            life.boundaryEnd = life.boundaryEnd.AddDays(-1);
        }
        else
        {
            life.boundaryStart = life.boundaryStart.AddDays(1);
            life.boundaryEnd = life.boundaryEnd.AddDays(1);
        }

        return life;
    }
    // 重载加法运算符
    public static Life operator +(Life life1, Life life2)
    {
        long newBoundaryStart = Mathf.Min(life1.boundaryStart, life2.boundaryStart);
        long newBoundaryEnd = Mathf.Max(life1.boundaryEnd, life2.boundaryEnd);

        Life newLife = new Life(life1.lifespan, newBoundaryStart, newBoundaryEnd);
        return newLife;
    }
}